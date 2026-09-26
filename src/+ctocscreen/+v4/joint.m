function [best,report,warm]=joint(candidate,ids,theta,scope,focus,eph,c,budget,resume)
%JOINT Restoration and fuel optimization with actual full-history replay.
if nargin<9, resume=[]; end
clock=tic; best=[]; warm=[];
report=struct('scope',scope,'status','not_started','iterations',0,'accepted_steps',0, ...
 'solver_seconds',0,'evaluation_seconds',0,'replay_seconds',0,'seconds',0, ...
 'solver_exitflags',[],'rebase_count',0,'stall_rebases',0,'replay_count',0,'initial_V',NaN,'final_V',NaN, ...
 'initial_J',candidate.actual.total_dv_km_s,'final_J',NaN,'opened_initial',false, ...
 'opened_burns',[],'control_changed',false,'active_passed',false,'resumed',false, ...
 'step_diagnostics',{{}},'last_replay_status','','last_replay_failure_id','','failure_id','','failure_reason','');
lastEvaluated=[]; lastReplayZ=[]; iterCost=0;
% Until a replay is timed here, estimate it from the measured ~0.05 s per arc
% of run 04 (fixed-control replay with 35-target scans).
replayCost=0; if isfield(candidate,'trace'), replayCost=.05*numel(candidate.trace.arcs); end
try
 p=ctocscreen.v4.problem(candidate,ids,theta,scope,focus,eph,c); z=p.z0;
 report.opened_initial=all(p.free(p.ix(:,1))); report.opened_burns=p.opened_burns;
 fuel=acceptable(candidate.actual,ids,c); if fuel, best=candidate; end
 trust=c.trust_initial; lambda=c.penalty_initial; stall=0; e=[];
 taskKey=[ctocscreen.v4.controlKey(candidate.q),scope,sprintf('%.17g,',[ids(:);theta(:);focus;p.N])];
 resumeCount=0;
 if ~isempty(resume)
  assert(strcmp(taskKey,resume.task_key)&&numel(resume.z)==p.n,'ctocscreen:v4:resumeTask','Recovery identity changed.');
  z=resume.z; trust=resume.trust; lambda=resume.penalty; fuel=resume.fuel;
  best=resume.best; report.resumed=true; resumeCount=resume.resume_count+1;
 end
 while report.iterations<c.joint_iterations
  % The allowance covers the next iteration and the closing replay, both
  % measured in this call; the first iteration always runs.
  if report.iterations>0&&toc(clock)+iterCost+replayCost>=budget, report.status='budget_interrupted'; break; end
  it0=tic;
  if isempty(e)
   et=tic; e=ctocscreen.v4.evaluate(p,z,true); report.evaluation_seconds=report.evaluation_seconds+toc(et);
  end
  if isnan(report.initial_V), report.initial_V=e.V; end
  report.final_V=e.V;
  % Fuel steps require the model to be near feasible as well (spec 6.3). A
  % replay-accepted incumbent may lie between the model and acceptance limits;
  % it stays the incumbent while the model is first restored inside its margin.
  if fuel&&~e.near_feasible, fuel=false; end
  step=ctocscreen.v4.coneStep(p,z,e,trust,lambda,fuel,max(.01,budget-toc(clock)-replayCost));
  report.iterations=report.iterations+1; report.solver_seconds=report.solver_seconds+step.seconds;
  report.solver_exitflags(end+1)=step.exitflag;
  sd=struct('trust',trust,'exitflag',step.exitflag,'ok',step.ok,'linear_feasibility',Inf,'ratio',NaN,'predicted_drop',NaN);
  if isfield(step,'feasibility'), sd.linear_feasibility=step.feasibility; end
  if ~step.ok
   report.step_diagnostics{end+1}=sd;
   trust=trust/2; stall=stall+1; report.status='cone_failed'; iterCost=toc(it0);
   if trust<c.trust_min||stall>=8, break; end
   continue
  end
  [prediction,~,~]=ctocscreen.v4.modelMerit(p,z,e,step.d,lambda,fuel);
  current=e.V; if fuel, current=e.J/p.sv+lambda*e.V; end
  pred=current-prediction;
  sd.predicted_drop=pred;
  if pred<=1e-9*max(1,abs(current))
   report.step_diagnostics{end+1}=sd;
   % The model cannot improve further. If it is near feasible but this iterate
   % was never replayed, let the actual trajectory decide (and rebase on failure).
   if ~fuel&&e.near_feasible&&~isequal(z,lastReplayZ)&&report.rebase_count<c.rebase_limit
    report.stall_rebases=report.stall_rebases+1; checkReplay(true); iterCost=toc(it0);
    continue
   end
   report.status='restoration_stalled'; if fuel, report.status='local_stationary'; end
   break
  end
  trial=z+step.d; et=tic;
  try
   te=ctocscreen.v4.evaluate(p,trial,false); nonlinear=te.V;
   if ~fuel, nonlinear=nonlinear+c.restoration_step_weight*norm(step.d); end
   if fuel, nonlinear=te.J/p.sv+lambda*te.V; end
   ratio=(current-nonlinear)/pred;
  catch err
   ratio=-Inf; report.failure_id=err.identifier; report.failure_reason=err.message;
  end
  report.evaluation_seconds=report.evaluation_seconds+toc(et);
  sd.ratio=ratio; report.step_diagnostics{end+1}=sd;
  if ratio<.1
   trust=trust/2; stall=stall+1;
  else
   z=trial; e=[]; report.accepted_steps=report.accepted_steps+1; stall=0; report.final_V=te.V;
   if ratio>.8&&norm(step.d,Inf)>.8*trust, trust=min(c.trust_max,trust*1.6); end
   if ratio<.25, trust=trust/1.5; end
   if fuel&&step.slack>c.restore_tolerance, lambda=min(c.penalty_max,lambda*5); end
   % Fuel iterates are nearly always inside the model gap, so they replay on the
   % three-step cadence; restoration replays as soon as the model is near feasible.
   if (~fuel&&te.near_feasible)||mod(report.accepted_steps,3)==0, checkReplay(te.near_feasible); end
  end
  iterCost=toc(it0);
  if trust<c.trust_min, report.status='trust_stalled'; break; end
 end
 if report.iterations>=c.joint_iterations&&~ismember(report.status,{'local_stationary','restoration_stalled','trust_stalled','cone_failed','budget_interrupted'})
  report.status='iteration_limit';
 end
 if isempty(lastEvaluated)||~isequal(z,lastReplayZ)
  rt=tic; q=ctocscreen.v4.decode(p,z); [a,tr]=ctocscreen.v4.replay(q,eph,c,false,candidate);
  report.replay_seconds=report.replay_seconds+toc(rt); report.replay_count=report.replay_count+1;
  report.last_replay_status=a.status; report.last_replay_failure_id=a.failure_id;
  trialCandidate=candidate; trialCandidate.q=q; trialCandidate.actual=a; trialCandidate.trace=tr;
  lastEvaluated=trialCandidate;
  if acceptable(a,ids,c)&&(isempty(best)||a.total_dv_km_s<best.actual.total_dv_km_s-1e-9), best=trialCandidate; end
 end
 if isempty(best)&&~isempty(lastEvaluated)&&~strcmp(lastEvaluated.actual.status,'propagation_failure')
  best=lastEvaluated;
 end
 if ~isempty(best)
  best.q=ctocscreen.v4.canonical(best.q,eph.model); best.q.witness=best.actual.witness_times_s;
  report.final_J=best.actual.total_dv_km_s;
  report.active_passed=acceptable(best.actual,ids,c);
  report.control_changed=~strcmp(ctocscreen.v4.controlKey(best.q),ctocscreen.v4.controlKey(candidate.q));
  % A changed trajectory needs its own remaining-cost estimate.
  if report.control_changed&&isfield(best,'heuristic_H'), best.heuristic_H=NaN; end
 end
 warm=struct('z',z,'task_ids',ids,'theta',theta,'scope',scope,'focus',focus,'trust',trust,'penalty',lambda, ...
  'physical_key',ctocscreen.v4.controlKey(candidate.q),'kind','restoration_iterate', ...
  'task_key',taskKey,'candidate',candidate,'best',best,'fuel',fuel,'resume_count',resumeCount, ...
  'resumable',report.rebase_count==0);
 if toc(clock)>=budget, report.status='budget_interrupted'; end
catch err
 report.status='numerical_failure'; report.failure_id=err.identifier; report.failure_reason=err.message;
end
report.seconds=toc(clock);
 function checkReplay(nearFeasible)
  % Fixed-control replay decides acceptance; a failed near-feasible model is
  % rebased on the actual states so hidden residuals become visible.
  rt=tic; q=ctocscreen.v4.decode(p,z); [a,tr]=ctocscreen.v4.replay(q,eph,c,false,candidate);
  spent=toc(rt); replayCost=spent;
  report.replay_seconds=report.replay_seconds+spent; report.replay_count=report.replay_count+1;
  report.last_replay_status=a.status; report.last_replay_failure_id=a.failure_id;
  trialCandidate=candidate; trialCandidate.q=q; trialCandidate.actual=a; trialCandidate.trace=tr;
  lastEvaluated=trialCandidate; lastReplayZ=z;
  if acceptable(a,ids,c)
   if isempty(best)||a.total_dv_km_s<best.actual.total_dv_km_s-1e-9, best=trialCandidate; end
   fuel=true; report.status='actual_feasible';
  elseif nearFeasible&&~strcmp(a.status,'propagation_failure')&&report.rebase_count<c.rebase_limit
   % Rebase auxiliary nodes only; the physical controls remain untouched.
   p2=ctocscreen.v4.problem(trialCandidate,ids,q.witness(ids),scope,focus,eph,c);
   if ~a.height_passed, p2.height_fractions=(1:7)/8; end
   p=p2; z=p.z0; e=[]; fuel=false; report.rebase_count=report.rebase_count+1;
   lastReplayZ=z; report.status='rebased_on_replay';
  else
   % Restore nonlinear connection defects before taking more fuel steps.
   fuel=false; report.status='restoring_actual_coverage';
  end
 end
end
function yes=acceptable(a,ids,c)
yes=~strcmp(a.status,'propagation_failure')&&a.initial_passed&&a.height_passed ...
 &&all(a.distance_km(ids)<=c.search_radius_km);
end
