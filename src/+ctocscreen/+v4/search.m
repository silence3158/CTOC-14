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
checkpoint=[]; checkpointElapsed=NaN; checkpointReport=[]; absorbTabu={}; tailDone={};
nextRoot=1; rootServices=zeros(1,16); stats=struct('rounds',0,'root_count',0,'expanded',0,'generated',0, ...
 'shared_calls',0,'full_calls',0,'structure_calls',0,'resume_calls',0,'joint_iterations',0, ...
 'actual_joint_improvements',0,'suffix_rebuilds',0,'independent_checks',0,'complete_found',0, ...
 'first_complete_s',NaN,'first_complete_dv_km_s',NaN,'notification',false, ...
 'expansion_seconds',0,'joint_seconds',0,'estimate_calls',0,'estimate_seconds',0,'absorb_calls',0,'absorbed',0,'tail_calls',0,'b_seconds',0,'verification_seconds',0,'source_unchanged',false,'stopped_on_complete',false);
try
for k=1:c.root_count
 node=ctocscreen.v4.root(nextRoot,eph,c,stream); nextRoot=nextRoot+1;
 beam{end+1}=node; stats.root_count=stats.root_count+1;
end
deadline=max(0,c.budget_s-c.verify_reserve_s); best=beam{1};
record('start',struct('cold_start',true,'budget_s',c.budget_s,'seed',c.seed));
% 2026-09-26 schedule (docs/V4_A_LAYER_CROSSING_20260926.md, section 10):
% construction first; periodic full-history B and structure operators are off
% (E3: early-control steps are amplified 1e4-1e9 on long chains). Tail B and
% encounter absorption act only on the last burns.
bSeconds=0;
% stop_on_complete ends the search at its first independently verified
% complete solution (used as phase 1 of a separate experiment entry).
while toc(clock)<deadline&&~stopNow()
 stats.rounds=stats.rounds+1; iteration=stats.rounds;
 if isempty(beam)||mod(iteration,c.root_every)==0
  node=ctocscreen.v4.root(nextRoot,eph,c,stream); nextRoot=nextRoot+1;
  beam{end+1}=node; stats.root_count=stats.root_count+1; parentIndex=numel(beam);
 else
  % Rotate through the beam: least-expanded first, ties by estimated J+H.
  attempts=cellfun(@(n)n.attempts,beam);
  costs=cellfun(@(n)n.actual.total_dv_km_s+estimateOf(n),beam);
  open=cellfun(@(n)n.actual.visit_count<35&&n.q.T<eph.model.horizon_s-1,beam);
  if ~any(open), open=true(size(beam)); end
  key=attempts+1e-3*costs/max(1,max(costs)); key(~open)=Inf;
  [~,parentIndex]=min(key);
  if rand(stream)<c.exploration, pool=find(open); parentIndex=pool(randi(stream,numel(pool))); end
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
 if stopNow(), break; end
 % Step 3: fold a following encounter into the newest burn (one burn, two targets).
 transfer=cellfun(@(n)any(strcmp(n.origin,{'crossing_transfer','delayed_impulse','absorbed_encounter'})),children);
 for j=find(transfer)
  if toc(clock)>=deadline, break; end
  bt=tic; [child,rr]=ctocscreen.v4.absorb(children{j},eph,c,min(c.absorb_seconds,deadline-toc(clock)),absorbTabu);
  bSeconds=bSeconds+toc(bt); stats.absorb_calls=stats.absorb_calls+1;
  absorbTabu{end+1}=sprintf('%s|%d',ctocscreen.v4.controlKey(children{j}.q),rr.target); %#ok<AGROW>
  if ~isempty(child), stats.absorbed=stats.absorbed+1; children{end+1}=child; consider(child); end %#ok<AGROW>
  record('absorb',rmfield(rr,'joint'));
 end
 % Step 2: tail B on the cheapest deep candidate while B stays under its share.
 if bSeconds<c.b_share*toc(clock)&&toc(clock)<deadline
  pool=[beam,children]; deep=cellfun(@(n)numel(n.q.tau)>=2&&n.actual.visit_count>=2,pool);
  if any(deep)
   pool=pool(deep); counts=cellfun(@(n)n.actual.visit_count,pool);
   top=pool(counts==max(counts)); [~,j]=min(cellfun(@(n)n.actual.total_dv_km_s,top)); seed=top{j};
   key=[ctocscreen.v4.controlKey(seed.q),'|tail'];
   if ~any(strcmp(tailDone,key))
    tailDone{end+1}=key; ids=find(seed.actual.distance_km<=1); M=numel(seed.q.tau); %#ok<AGROW>
    cc=c; cc.joint_iterations=c.tail_iterations; bt=tic;
    [child,jr,wr]=ctocscreen.v4.joint(seed,ids,seed.actual.witness_times_s(ids),'tail', ...
     seed.q.tau(max(1,M-c.tail_burns+1)),eph,cc,min(c.tail_seconds,deadline-toc(clock)));
    bSeconds=bSeconds+toc(bt); stats.tail_calls=stats.tail_calls+1; addJoint(jr,[]);
    if ~isempty(child)&&jr.active_passed&&jr.control_changed
     child.origin='tail_b'; child.attempts=0; children{end+1}=child; consider(child);
    end
   end
  end
 end
 children=children(cellfun(@(n)~strcmp(n.actual.status,'propagation_failure') ...
  &&n.actual.initial_passed&&n.actual.height_passed&&isfinite(n.actual.total_dv_km_s),children));
 et=tic; [beam,estimated]=ctocscreen.v4.selectBeam([beam,children],c,eph,stream);
 stats.estimate_calls=stats.estimate_calls+estimated; stats.estimate_seconds=stats.estimate_seconds+toc(et);
 stats.b_seconds=bSeconds;
 if mod(iteration,c.log_every)==0
  fprintf('V4 round=%d elapsed=%.1f best=%d/35 dv=%.9f beam=%d absorb=%d/%d tail=%d B=%.0fs\n', ...
   iteration,toc(clock),best.actual.visit_count,best.actual.total_dv_km_s,numel(beam), ...
   stats.absorbed,stats.absorb_calls,stats.tail_calls,bSeconds);
 end
 record('round',struct('round',iteration,'best_visits',best.actual.visit_count, ...
  'best_dv',best.actual.total_dv_km_s,'beam_size',numel(beam)));
end
stats.search_seconds=toc(clock); stats.stopped_on_complete=stopNow();
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
 function yes=stopNow()
  yes=c.stop_on_complete>0&&~isempty(bestVerified);
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
