function [best,report,warm]=joint(candidate,ids,theta,scope,focus,eph,c,budget,resume)
%JOINT Restoration and fuel optimization with actual full-history replay.
if nargin<9, resume=[]; end
clock=tic; best=[]; warm=[];
report=struct('scope',scope,'status','not_started','iterations',0,'accepted_steps',0, ...
 'solver_seconds',0,'evaluation_seconds',0,'replay_seconds',0,'seconds',0, ...
 'solver_exitflags',[],'rebase_count',0,'initial_V',NaN,'final_V',NaN, ...
 'initial_J',candidate.actual.total_dv_km_s,'final_J',NaN,'opened_initial',false, ...
 'opened_burns',[],'control_changed',false,'active_passed',false,'resumed',false,'failure_id','','failure_reason','');
try
 p=ctocscreen.v4.problem(candidate,ids,theta,scope,focus,eph,c); z=p.z0;
 report.opened_initial=all(p.free(p.ix(:,1))); report.opened_burns=p.opened_burns;
  fuel=acceptable(candidate.actual,ids,c); if fuel, best=candidate; end
 trust=c.trust_initial; lambda=c.penalty_initial; stall=0; lastEvaluated=[]; lastReplayZ=[];
 taskKey=[ctocscreen.v4.controlKey(candidate.q),scope,sprintf('%.17g,',[ids(:);theta(:);focus;p.N])];
 resumeCount=0;
 if ~isempty(resume)
  assert(strcmp(taskKey,resume.task_key)&&numel(resume.z)==p.n,'ctocscreen:v4:resumeTask','Recovery identity changed.');
  z=resume.z; trust=resume.trust; lambda=resume.penalty; fuel=resume.fuel;
  best=resume.best; report.resumed=true; resumeCount=resume.resume_count+1;
 end
 while report.iterations<c.joint_iterations&&toc(clock)<budget
  et=tic; e=ctocscreen.v4.evaluate(p,z,true); report.evaluation_seconds=report.evaluation_seconds+toc(et);
  if isnan(report.initial_V), report.initial_V=e.V; end
  report.final_V=e.V;
  if toc(clock)>=budget, report.status='budget_interrupted'; break; end
  step=ctocscreen.v4.coneStep(p,z,e,trust,lambda,fuel,budget-toc(clock));
  report.iterations=report.iterations+1; report.solver_seconds=report.solver_seconds+step.seconds;
  report.solver_exitflags(end+1)=step.exitflag;
  if ~step.ok
   trust=trust/2; stall=stall+1; report.status='cone_failed';
   if trust<c.trust_min||stall>=3, break; end
   continue
  end
  [prediction,~,~]=ctocscreen.v4.modelMerit(p,z,e,step.d,lambda,fuel);
  current=e.V; if fuel, current=e.J/p.sv+lambda*e.V; end
  pred=current-prediction;
  if pred<=1e-9*max(1,abs(current))
   report.status='restoration_stalled'; if fuel, report.status='local_stationary'; end
   break
  end
  trial=z+step.d; et=tic;
  try
   te=ctocscreen.v4.evaluate(p,trial,false); nonlinear=te.V;
   if fuel, nonlinear=te.J/p.sv+lambda*te.V; end
   ratio=(current-nonlinear)/pred;
  catch err
   ratio=-Inf; report.failure_id=err.identifier; report.failure_reason=err.message;
  end
  report.evaluation_seconds=report.evaluation_seconds+toc(et);
  if ratio<.1
   trust=trust/2; stall=stall+1;
  else
   z=trial; report.accepted_steps=report.accepted_steps+1; stall=0; report.final_V=te.V;
   if ratio>.8&&norm(step.d,Inf)>.8*trust, trust=min(c.trust_max,trust*1.6); end
   if ratio<.25, trust=trust/1.5; end
   if fuel&&step.slack>c.restore_tolerance, lambda=min(c.penalty_max,lambda*5); end
   if te.near_feasible||mod(report.accepted_steps,3)==0
    rt=tic; q=ctocscreen.v4.decode(p,z); [a,tr]=ctocscreen.v4.replay(q,eph,c);
    report.replay_seconds=report.replay_seconds+toc(rt);
    trialCandidate=candidate; trialCandidate.q=q; trialCandidate.actual=a; trialCandidate.trace=tr;
    lastEvaluated=trialCandidate;
    lastReplayZ=z;
    if acceptable(a,ids,c)
     if isempty(best)||a.total_dv_km_s<best.actual.total_dv_km_s-1e-9, best=trialCandidate; end
     fuel=true; report.status='actual_feasible';
    elseif te.near_feasible&&~strcmp(a.status,'propagation_failure')&&report.rebase_count<2
     % Rebase auxiliary nodes only; the physical controls remain untouched.
     p2=ctocscreen.v4.problem(trialCandidate,ids,q.witness(ids),scope,focus,eph,c);
     if ~a.height_passed, p2.height_fractions=(1:7)/8; end
     p=p2; z=p.z0; fuel=false; report.rebase_count=report.rebase_count+1;
    end
   end
  end
  if trust<c.trust_min||stall>=4, report.status='trust_stalled'; break; end
 end
 if isempty(lastEvaluated)||~isequal(z,lastReplayZ)
  rt=tic; q=ctocscreen.v4.decode(p,z); [a,tr]=ctocscreen.v4.replay(q,eph,c);
  report.replay_seconds=report.replay_seconds+toc(rt);
  trialCandidate=candidate; trialCandidate.q=q; trialCandidate.actual=a; trialCandidate.trace=tr;
  lastEvaluated=trialCandidate;
  if acceptable(a,ids,c)&&(isempty(best)||a.total_dv_km_s<best.actual.total_dv_km_s-1e-9), best=trialCandidate; end
 end
 if isempty(best)&&~isempty(lastEvaluated), best=lastEvaluated; end
 if ~isempty(best)
  best.q=ctocscreen.v4.canonical(best.q,eph.model); best.q.witness=best.actual.witness_times_s;
  report.final_J=best.actual.total_dv_km_s;
  report.active_passed=acceptable(best.actual,ids,c);
  report.control_changed=~strcmp(ctocscreen.v4.controlKey(best.q),ctocscreen.v4.controlKey(candidate.q));
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
end
function yes=acceptable(a,ids,c)
yes=~strcmp(a.status,'propagation_failure')&&a.initial_passed&&a.height_passed ...
 &&all(a.distance_km(ids)<=c.search_radius_km);
end
