function result=search(folder,eph,c)
%SEARCH Independent V4 A+B scheduling, cold roots, and fixed-control validation.
assert(~isfolder(folder),'ctocscreen:v4:existingRun','A new output folder is required.');
mkdir(folder); signature=ctocscreen.v4.signature(); clock=tic;
stream=RandStream('mt19937ar','Seed',c.seed); memory=[]; beam={};
manifest=struct('config',c,'signature',signature,'target_signature',eph.signature, ...
 'matlab_version',version,'cold_start',true,'initial_candidates',0,'resume_file','', ...
 'history_inputs',{{}},'started_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss')));
save(fullfile(folder,'manifest.mat'),'manifest');
fid=fopen(fullfile(folder,'events.jsonl'),'w','n','UTF-8'); assert(fid>0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
events={}; joints={}; warm={}; pending={}; best=[]; bestComplete=[]; bestVerified=[]; verifiedKeys={};
checkpoint=[]; checkpointElapsed=NaN; checkpointReport=[]; sharedTabu=struct('keys',{{}},'failures',[]); structureCount=0;
nextRoot=1; rootServices=zeros(1,16); stats=struct('rounds',0,'root_count',0,'expanded',0,'generated',0, ...
 'shared_calls',0,'full_calls',0,'structure_calls',0,'resume_calls',0,'joint_iterations',0, ...
 'actual_joint_improvements',0,'suffix_rebuilds',0,'independent_checks',0,'complete_found',0, ...
 'first_complete_s',NaN,'first_complete_dv_km_s',NaN,'notification',false, ...
 'expansion_seconds',0,'joint_seconds',0,'estimate_calls',0,'estimate_seconds',0,'verification_seconds',0,'source_unchanged',false);
try
for k=1:c.root_count
 node=ctocscreen.v4.root(nextRoot,eph,c,stream); nextRoot=nextRoot+1;
 beam{end+1}=node; stats.root_count=stats.root_count+1;
end
deadline=max(0,c.budget_s-c.verify_reserve_s); best=beam{1};
record('start',struct('cold_start',true,'budget_s',c.budget_s,'seed',c.seed));
while toc(clock)<deadline
 stats.rounds=stats.rounds+1; iteration=stats.rounds;
 if isempty(beam)||mod(iteration,c.root_every)==0
  node=ctocscreen.v4.root(nextRoot,eph,c,stream); nextRoot=nextRoot+1;
  beam{end+1}=node; stats.root_count=stats.root_count+1; parentIndex=numel(beam);
 else
  visits=cellfun(@(n)n.actual.visit_count,beam); attempts=cellfun(@(n)n.attempts,beam);
  % Estimated mission cost J+H (km/s) trades against depth; H is heuristic.
  costs=cellfun(@(n)n.actual.total_dv_km_s+estimateOf(n),beam);
  weights=ones(size(beam));
  for k=1:numel(beam), [memory,weights(k)]=ctocscreen.v4.feedback(memory,'query',beam{k},c); end
  score=visits-.8*attempts-c.parent_cost_weight*costs+log(weights);
  [~,parentIndex]=max(score);
  if rand(stream)<c.exploration
   probabilities=exp(score-max(score)); probabilities=probabilities/sum(probabilities);
   parentIndex=find(rand(stream)<=cumsum(probabilities),1);
  end
  if mod(iteration,3)==0
   roots=cellfun(@(n)n.root_id,beam); rootServices(max(roots))=rootServicesAt(max(roots));
   least=min(rootServices(roots)); eligible=find(rootServices(roots)==least);
   [~,j]=max(score(eligible)); parentIndex=eligible(j);
  end
 end
 parent=beam{parentIndex}; parent.attempts=parent.attempts+1; beam{parentIndex}=parent;
 rootServices(parent.root_id)=rootServicesAt(parent.root_id)+1;
 if parent.q.T>=eph.model.horizon_s-1
  [rebuilt,ar]=ctocscreen.v4.rebuild(parent,eph,c,parent.attempts);
  children={};
  if ~isempty(rebuilt)
   rebuilt.root_id=nextRoot; nextRoot=nextRoot+1; children={rebuilt}; stats.suffix_rebuilds=stats.suffix_rebuilds+1;
  end
  record('suffix_rebuild',ar);
 else
  [children,ar,memory]=ctocscreen.v4.expand(parent,eph,c,stream,memory,min(c.action_seconds,deadline-toc(clock)));
  stats.expanded=stats.expanded+1; stats.expansion_seconds=stats.expansion_seconds+ar.seconds;
  record('expand',struct('root_id',parent.root_id,'parent_visits',parent.actual.visit_count,'report',ar));
 end
 stats.generated=stats.generated+numel(children);
 for j=1:numel(children), consider(children{j}); end
 if ~isempty(pending)&&mod(iteration,3)==0&&toc(clock)<deadline
  job=pending{1}; pending(1)=[];
  [child,jr,wr]=ctocscreen.v4.joint(job.candidate,job.task_ids,job.theta,job.scope,job.focus,eph,c, ...
   min(c.scope_seconds(2),deadline-toc(clock)),job);
  stats.resume_calls=stats.resume_calls+1; addJoint(jr,wr);
  if ~isempty(child), children{end+1}=child; consider(child); end
 end
 % One shared-arc job; its scope escalates with executed failures of that task.
 if ~isempty(children)&&toc(clock)<deadline
  counts=cellfun(@(n)n.actual.visit_count,children); [~,j]=max(counts); selected=children{j};
  [child,jr,wr,sharedTabu]=ctocscreen.v4.shared(selected,eph,c,min(c.scope_seconds,deadline-toc(clock)),sharedTabu);
  stats.shared_calls=stats.shared_calls+1; addJoint(jr,wr);
  if ~isempty(child), children{end+1}=child; consider(child); end
 end
 if mod(iteration,c.full_every)==0&&toc(clock)<deadline
  representatives=[beam,children]; counts=cellfun(@(n)n.actual.visit_count,representatives);
  [~,j]=max(counts); seed=representatives{j}; ids=find(seed.actual.distance_km<=1);
  if ~isempty(ids)
   continuation=[]; key=ctocscreen.v4.controlKey(seed.q);
   for j=numel(pending):-1:1
    if strcmp(pending{j}.scope,'full')&&strcmp(pending{j}.physical_key,key) ...
      &&isequal(pending{j}.task_ids,ids)&&isequal(pending{j}.theta,seed.actual.witness_times_s(ids))
     continuation=pending{j}; pending(j)=[]; break
    end
   end
   [child,jr,wr]=ctocscreen.v4.joint(seed,ids,seed.actual.witness_times_s(ids),'full',seed.q.T,eph,c, ...
    min(c.scope_seconds(3),deadline-toc(clock)),continuation);
   stats.full_calls=stats.full_calls+1; addJoint(jr,wr);
   if ~isempty(child), child.origin='full_history'; child.heuristic_H=NaN; children{end+1}=child; consider(child); end
  end
 end
 if mod(iteration,c.structure_every)==0&&toc(clock)<deadline
  % A fresh root has no control history; use the best-covered nonempty trajectory.
  target=parent;
  if target.q.T<=0||target.q.T>=eph.model.horizon_s-1
   pool=[beam,children]; ok=cellfun(@(n)n.q.T>0&&n.q.T<eph.model.horizon_s-1,pool);
   target=[];
   if any(ok)
    pool=pool(ok); counts=cellfun(@(n)n.actual.visit_count,pool); [~,j]=max(counts); target=pool{j};
   end
  end
  if ~isempty(target)
   structureCount=structureCount+1;
   [child,jr,wr]=ctocscreen.v4.restructure(target,eph,c,structureCount,min(c.scope_seconds(2),deadline-toc(clock)));
   stats.structure_calls=stats.structure_calls+1; jr.structure_index=structureCount; addJoint(jr,wr);
   if ~isempty(child)
    if strcmp(child.origin,'suffix_rebuild')
     child.root_id=nextRoot; nextRoot=nextRoot+1; stats.suffix_rebuilds=stats.suffix_rebuilds+1;
    end
    children{end+1}=child; consider(child);
   end
  end
 end
 children=children(cellfun(@(n)~strcmp(n.actual.status,'propagation_failure') ...
  &&n.actual.initial_passed&&n.actual.height_passed&&isfinite(n.actual.total_dv_km_s),children));
 % Feedback affects competing physical children, with a nonzero exploration floor.
 for j=1:numel(children)
  [memory,weight]=ctocscreen.v4.feedback(memory,'query',children{j},c);
  if weight<1&&rand(stream)>max(c.exploration,weight)
   children{j}.attempts=children{j}.attempts+1;
  end
 end
 et=tic; [beam,estimated]=ctocscreen.v4.selectBeam([beam,children],c,eph);
 stats.estimate_calls=stats.estimate_calls+estimated; stats.estimate_seconds=stats.estimate_seconds+toc(et);
 [memory,~]=ctocscreen.v4.feedback(memory,'evaporate',[],c);
 if mod(iteration,c.log_every)==0
  fprintf('V4 round=%d elapsed=%.1f best=%d/35 dv=%.9f beam=%d B=%d/%d/%d\n', ...
   iteration,toc(clock),best.actual.visit_count,best.actual.total_dv_km_s,numel(beam), ...
   stats.shared_calls,stats.full_calls,stats.structure_calls);
 end
 record('round',struct('round',iteration,'best_visits',best.actual.visit_count, ...
  'best_dv',best.actual.total_dv_km_s,'beam_size',numel(beam)));
end
stats.search_seconds=toc(clock);
chosen=best; if ~isempty(bestComplete), chosen=bestComplete; end
if ~isempty(bestVerified), chosen=bestVerified; end
if ~isempty(bestVerified)&&strcmp(ctocscreen.v4.controlKey(chosen.q),ctocscreen.v4.controlKey(bestVerified.q))
 verification=bestVerified.actual; vtrace=bestVerified.trace;
else
 vt=tic; [verification,vtrace]=ctocscreen.v4.replay(chosen.q,eph,c,true);
 stats.verification_seconds=stats.verification_seconds+toc(vt); stats.independent_checks=stats.independent_checks+1;
end
if verification.passed
 verifiedNode=chosen; verifiedNode.actual=verification; verifiedNode.trace=vtrace;
 [memory,~]=ctocscreen.v4.feedback(memory,'observe',verifiedNode,c);
 if isempty(bestVerified)||verification.total_dv_km_s<bestVerified.actual.total_dv_km_s, bestVerified=verifiedNode; end
 if verification.total_dv_km_s<c.notify_dv_km_s, stats.notification=true; end
end
if ~isempty(checkpoint)&&c.budget_s>c.checkpoint_s
 vt=tic;
 if strcmp(ctocscreen.v4.controlKey(checkpoint.q),ctocscreen.v4.controlKey(chosen.q))
  checkpointReport=verification;
 else
  [checkpointReport,~]=ctocscreen.v4.replay(checkpoint.q,eph,c,true);
  stats.independent_checks=stats.independent_checks+1;
 end
 stats.verification_seconds=stats.verification_seconds+toc(vt);
 snapshot=struct('candidate',thin(checkpoint),'verification',checkpointReport, ...
  'candidate_elapsed_s',checkpointElapsed,'limit_s',c.checkpoint_s,'verified_elapsed_s',toc(clock));
 save(fullfile(folder,'budget_checkpoint.mat'),'snapshot');
 record('budget_checkpoint',struct('limit_s',c.checkpoint_s,'candidate_elapsed_s',checkpointElapsed, ...
  'visits',checkpointReport.visit_count,'dv',checkpointReport.total_dv_km_s,'passed',checkpointReport.passed));
end
stats.source_unchanged=isequaln(signature,ctocscreen.v4.signature());
assert(stats.source_unchanged,'ctocscreen:v4:changedSource','Sources changed during this experiment.');
record('final_verification',verification);
result=struct('manifest',manifest,'stats',stats,'best',thin(chosen),'verification',verification, ...
 'best_complete',thin(bestComplete),'best_verified',thin(bestVerified),'feedback',memory, ...
 'events',{events},'joint_reports',{joints},'recovery_states',{warm},'rng_state',stream.State, ...
 'checkpoint_verification',checkpointReport,'checkpoint_candidate_elapsed_s',checkpointElapsed);
result.stats.total_seconds=toc(clock);
save(fullfile(folder,'result.mat'),'result','-v7.3');
if ~isempty(bestVerified)
 verified=thin(bestVerified); save(fullfile(folder,'verified_complete.mat'),'verified');
 if verified.actual.total_dv_km_s<=c.search_max_dv_km_s
  save(fullfile(folder,'elite.mat'),'verified');
 else
  save(fullfile(folder,'rejected_complete.mat'),'verified');
 end
end
ctocscreen.v4.writeReport(folder,result,vtrace,eph);
result.stats.total_seconds=toc(clock); save(fullfile(folder,'result.mat'),'result','-v7.3');
fprintf('V4 FINAL independent=%d/35 dv=%.12f passed=%d elapsed=%.3f folder=%s\n', ...
 verification.visit_count,verification.total_dv_km_s,verification.passed,result.stats.total_seconds,folder);
catch err
 interrupted=struct('status','runtime_interrupted','failure_id',err.identifier,'failure_reason',err.message, ...
  'elapsed_s',toc(clock),'best',thin(best),'best_verified',thin(bestVerified),'manifest',manifest,'stats',stats);
 save(fullfile(folder,'interrupted.mat'),'interrupted','-v7.3');
 rethrow(err);
end
 function record(event,payload)
  row=struct('event',event,'elapsed_s',toc(clock),'payload',payload);
  events{end+1}=row; fprintf(fid,'%s\n',jsonencode(row));
 end
 function addJoint(jr,wr)
  joints{end+1}=jr;
  if ~isempty(wr)
   warm{end+1}=rmfield(wr,{'candidate','best'});
   if wr.resumable&&wr.resume_count<8&&~ismember(jr.status,{'local_stationary','restoration_stalled','numerical_failure','trust_stalled'})
    duplicate=cellfun(@(j)strcmp(j.task_key,wr.task_key),pending); pending(duplicate)=[];
    pending{end+1}=wr; if numel(pending)>12, pending(1)=[]; end
   end
  end
  if isfield(jr,'iterations'), stats.joint_iterations=stats.joint_iterations+jr.iterations; end
  stats.joint_seconds=stats.joint_seconds+jr.seconds;
  if isfield(jr,'final_J')&&jr.final_J<jr.initial_J-1e-7&&jr.control_changed&&jr.active_passed
   stats.actual_joint_improvements=stats.actual_joint_improvements+1;
  end
  record('joint',jr);
 end
 function value=estimateOf(n)
  value=0; if isfield(n,'heuristic_H')&&~isnan(n.heuristic_H), value=n.heuristic_H; end
 end
 function value=rootServicesAt(id)
  value=0; if id<=numel(rootServices), value=rootServices(id); end
 end
 function consider(node)
  if ~node.actual.initial_passed||~node.actual.height_passed||strcmp(node.actual.status,'propagation_failure'), return; end
  [memory,~]=ctocscreen.v4.feedback(memory,'observe',node,c);
  if isempty(best)||node.actual.visit_count>best.actual.visit_count ...
    ||(node.actual.visit_count==best.actual.visit_count&&node.actual.total_dv_km_s<best.actual.total_dv_km_s)
   best=node; record('best',struct('visits',node.actual.visit_count,'dv',node.actual.total_dv_km_s, ...
    'root_id',node.root_id,'origin',node.origin));
   if c.checkpoint_s>0&&toc(clock)<=c.checkpoint_s
    checkpoint=best; checkpointElapsed=toc(clock);
   end
  end
  if ~node.actual.passed, return; end
  if isempty(bestComplete)||node.actual.total_dv_km_s<bestComplete.actual.total_dv_km_s, bestComplete=node; end
  key=ctocscreen.v4.controlKey(node.q);
  if any(strcmp(verifiedKeys,key)), return; end
  verifiedKeys{end+1}=key; vt0=tic;
  [verified,trv]=ctocscreen.v4.replay(node.q,eph,c,true);
  stats.verification_seconds=stats.verification_seconds+toc(vt0); stats.independent_checks=stats.independent_checks+1;
  if verified.passed
   stats.complete_found=stats.complete_found+1;
   if isnan(stats.first_complete_s)
    stats.first_complete_s=toc(clock); stats.first_complete_dv_km_s=verified.total_dv_km_s;
   end
   node.actual=verified; node.trace=trv;
   if isempty(bestVerified)||verified.total_dv_km_s<bestVerified.actual.total_dv_km_s, bestVerified=node; end
   [memory,~]=ctocscreen.v4.feedback(memory,'observe',node,c);
   if verified.total_dv_km_s<c.notify_dv_km_s
    stats.notification=true;
    fprintf('V4 USER_ACCEPTANCE_READY independent=35/35 raw_dv=%.12f time=%.3f\n',verified.total_dv_km_s,toc(clock));
    witness=thin(node); save(fullfile(folder,'acceptance_ready.mat'),'witness','manifest');
   end
  end
  record('independent',verified);
 end
end
function s=thin(node)
s=[]; if isempty(node), return; end
s=rmfield(node,'trace');
end
