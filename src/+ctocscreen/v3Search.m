function state=v3Search(folder,eph,c)
%V3SEARCH Beam ACO + unified shooting + feasible SA + independent export gate.
assert(~isfolder(folder),'ctocscreen:v3:folder','Use a new output folder.'); mkdir(folder);
% Verify actual solver license; no silent substitution by another algorithm.
probe=optimoptions('fmincon','Display','off');
fmincon(@(x)x*x,1,[],[],[],[],0,2,[],probe);
lsqlin(eye(2),[1;2],[],[],[],[],[-3;-3],[3;3],[],optimoptions('lsqlin','Display','off'));
if strcmp(c.fuel_algorithm,'sparse_sqp')
 quadprog(speye(2),[-1;-1],[],[],[],[],[-2;-2],[2;2],[],optimoptions('quadprog','Display','off'));
end
signature=ctocscreen.v3ImplementationSignature(eph); stream=RandStream('mt19937ar','Seed',c.seed);
for initial=1:numel(c.initial_candidates)
 c.initial_candidates{initial}=ctocscreen.v3Normalize(c.initial_candidates{initial},eph.model);
end
state=struct('config',c,'signature',signature,'iteration',0,'history',{{}},'failures',{{}}, ...
 'elite',[],'work_pool',{{}},'recovery_pool',{{}},'root_archive',{{}}, ...
 'temperature',c.temperature_km_s,'temperature_samples',[],'stalled',0, ...
 'elapsed_s',0,'budget_overrun_s',0,'export_archive',{{}}, ...
 'pending_initial',{c.initial_candidates},'dataset_kind','competition_targets_nominal_model', ...
 'construction_history',{{}},'best_partial',[],'pending_roots',{{}}, ...
 'lineage_keys',{{}},'lineage_visits',[],'lineage_stalls',[], ...
 'refined_keys',{{}},'duplicate_joint_skips',0,'reinforcement_ledger',{{}},'feedback_history',{{}}, ...
 'pool_history',{{}},'mutation_history',{{}},'fresh_queue',{{}},'replan_history',{{}}, ...
 'replan_queue',{{}},'replanned_keys',{{}},'beam_expansion_history',{{}}, ...
 'source_history',{{}},'colony_iterations',zeros(1,c.colony_count), ...
 'prefix_population',{{}},'prefix_feedback_history',{{}},'full_construction_history',{{}}, ...
 'protected_attempts',0,'prefix_restart_history',{{}}, ...
 'work_lineage_keys',{{}},'work_lineage_attempts',[],'policy_rejections',{{}},'policy_rejection_count',0, ...
 'policy_rejected_keys',{{}},'best_rejected_complete',[]);
if isfield(eph,'synthetic')&&eph.synthetic, state.dataset_kind='synthetic_test_only'; end
pheromone=containers.Map('KeyType','char','ValueType','double');
colonies=cell(1,c.colony_count);
for colony=1:c.colony_count
 colonies{colony}=ctocscreen.v3Roots(eph,c,stream);
 state.root_archive=[state.root_archive colonies{colony}];
 state.pending_roots=[state.pending_roots colonies{colony}];
end
active=1; beam=colonies{active};
% Initial inputs are handled once through pending_initial; do not also
% recycle verified complete inputs through the infeasible recovery queue.
if ~isempty(c.resume_file)
 old=load(c.resume_file,'state'); old=old.state;
 assert(isequaln(old.signature,signature),'ctocscreen:v3:staleCheckpoint','Checkpoint dependencies changed.');
 assert(old.config.seed==c.seed,'ctocscreen:v3:seed','Resume seed must match.');
 assert(old.config.colony_count==c.colony_count,'ctocscreen:v3:colonies','Resume colony count must match.');
 state=old; state.config=c; stream.State=state.rng_state; beam=state.beam; colonies=state.colonies;
 if ~isempty(state.pheromone_keys), pheromone=containers.Map(state.pheromone_keys,state.pheromone_values); end
 if ~isempty(state.elite)
  elite=state.elite; save(fullfile(folder,'elite.mat'),'elite','-v7.3');
 end
end
clock=tic; baseElapsed=state.elapsed_s;
state.initial_seed_checks={};
if isempty(c.resume_file)
 for initial=1:numel(state.pending_initial)
  if toc(clock)>=c.budget_s, break; end
  input=state.pending_initial{initial};
  if sum(vecnorm(input.delta_v_km_s,2,2))>c.search_max_dv_km_s, rejectSchedule(input,'initial'); continue; end
  if ~all(isfinite(input.witness_times_s)), continue; end
  checked=ctocscreen.v3Verify(input,eph,c); state.initial_seed_checks{end+1}=checked;
  if checked.passed
   input.witness_times_s=checked.witness_times_s;
   item=deposit(input,checked.total_dv_km_s);
   [state.work_pool,admission]=ctocscreen.v3SelectWorkPool(state.work_pool,item,c,false);
   state.pool_history{end+1}=admission;
   if isempty(state.elite)||checked.total_dv_km_s<state.elite.verification.total_dv_km_s
    state.elite=struct('schedule',input,'verification',checked,'elapsed_s',baseElapsed+toc(clock));
    state.export_archive{end+1}=state.elite;
    elite=state.elite; save(fullfile(folder,'elite.mat'),'elite','-v7.3');
   end
  end
 end
end
while toc(clock)<c.budget_s && (isempty(c.stop_file)||~isfile(c.stop_file))
 state.iteration=state.iteration+1; iter=state.iteration;
 active=mod(iter-1,c.colony_count)+1; beam=colonies{active};
 state.colony_iterations(active)=state.colony_iterations(active)+1;
 [~,pheromone]=ctocscreen.v3Pheromone('evaporate',pheromone,'',0,c);
 if mod(iter,c.expand_after)==0
  c.max_maneuvers=c.max_maneuvers+8; c.pulse_component_bound=c.pulse_component_bound*1.5;
  c.guided_max_revolutions=c.guided_max_revolutions+2;
  c.cost_lookahead_s=min(eph.model.horizon_s,c.cost_lookahead_s*1.5);
  c.action_durations_s=unique([c.action_durations_s,min(864000,max(c.action_durations_s)*1.5)]);
  fresh=ctocscreen.v3Roots(eph,c,stream); beam=ctocscreen.v3SelectBeam([beam fresh],c.beam_width,c);
  state.root_archive=[state.root_archive fresh];
  state.pending_roots=[state.pending_roots fresh];
 end
 protected=mod(state.colony_iterations(active)-1,c.protected_root_every)==0;
 if protected&&isempty(state.pending_roots)
  fresh=ctocscreen.v3Roots(eph,c,stream); state.root_archive=[state.root_archive fresh];
  state.pending_roots=fresh;
 end
 if protected
  state.protected_attempts=state.protected_attempts+1;
  if ~isempty(state.prefix_population)&&mod(state.protected_attempts,c.prefix_restart_every)==0
   member=state.prefix_population{randi(stream,numel(state.prefix_population))};
   node=ctocscreen.v3PrefixRestart(member,eph);
   state.root_archive{end+1}=node;
   state.prefix_restart_history{end+1}=struct('iteration',iter,'member_key',member.physical_key, ...
    'root_key',member.root,'source','current_run_prefix_population');
  else
   node=state.pending_roots{1}; state.pending_roots(1)=[];
  end
  if isfield(node,'seed_wait_s')&&node.seed_wait_s>0&&node.t==0
   try
    node=ctocscreen.v3ApplyAction(node,0,zeros(3,1),node.seed_wait_s,0,eph,c);
   catch err
    state.failures{end+1}=struct('id',err.identifier,'message',err.message,'stage','initial_coast');
   end
  end
  completionInputs={node};
 else
 expansionClock=tic; expansionConfig=c; expansionConfig.budget_s=min(c.beam_seconds,c.budget_s-toc(clock));
 for depth=1:c.beam_expansions
  nodes={};
  for b=1:numel(beam)
   branchConfig=expansionConfig;
   branchConfig.budget_s=toc(expansionClock)+max(0,expansionConfig.budget_s-toc(expansionClock))/(numel(beam)-b+1);
   [child,failed,expansion]=ctocscreen.v3Expand(beam{b},eph,branchConfig,stream,pheromone,expansionClock);
   expansion.iteration=iter; expansion.colony=active; expansion.depth=depth;
   state.beam_expansion_history{end+1}=expansion;
   nodes=[nodes child]; state.failures=[state.failures failed]; %#ok<AGROW>
   if toc(expansionClock)>=expansionConfig.budget_s, break; end
  end
  if ~isempty(nodes), beam=ctocscreen.v3SelectBeam(nodes,c.beam_width,c); end
  if toc(expansionClock)>=expansionConfig.budget_s, break; end
 end
 if toc(clock)>=c.budget_s, break; end
 % All selected prefixes continue through beam completion.
 completionInputs=beam;
 end
 cc=c; cc.completion_seconds=min(c.completion_seconds,max(.01,c.budget_s-toc(clock)));
 [node,construction,failed,constructed]=ctocscreen.v3Complete(completionInputs,eph,cc,stream,pheromone,true, ...
  c.cold_completion_enabled);
 construction.pruned_extensions=0;
 if construction.consistency.passed && construction.consistency.truncated
  [beam,construction.pruned_extensions]=ctocscreen.v3PruneUnstableBeam(beam,node);
 end
 construction.iteration=iter; construction.started_at_s=toc(clock)-construction.elapsed_s;
 construction.colony=active; construction.colony_iteration=state.colony_iterations(active);
 construction.protected_root=protected; construction.root_key=node.schedule.root_key;
 construction.root_method=node.root_method;
 state.construction_history{end+1}=construction; state.failures=[state.failures failed];
 beam=ctocscreen.v3SelectBeam([constructed,beam],c.beam_width,c);
 items={};
 for candidateIndex=1:numel(constructed)
 n=constructed{candidateIndex}; ids=find(n.visited).';
  if n.J>c.search_max_dv_km_s&&numel(ids)<35
   isNew=rejectSchedule(n.schedule,'construction_output',n.completion_check);
   if isNew, enqueue(n.schedule,1:35,[],c.recovery_stall_limit); end
   continue;
  end
  if numel(ids)<2||any(strcmp(state.refined_keys,ctocscreen.v3ScheduleKey(n.schedule))), continue; end
  if any(strcmp(state.policy_rejected_keys,ctocscreen.v3ScheduleKey(n.schedule))), continue; end
  items{end+1}=struct('schedule',n.schedule,'target_ids',ids,'node',n,'operator','beam_completion');
 end
 [state.fresh_queue,construction.queued_candidates]=ctocscreen.v3QueueCandidates(state.fresh_queue,items,c.fresh_queue_size);
 % A completed construction must not wait behind, or be dropped for, old prefixes.
 for j=1:numel(items)
  item=items{j};
  if numel(item.target_ids)==35&&~any(cellfun(@(p)strcmp(ctocscreen.v3ScheduleKey(p.schedule), ...
      ctocscreen.v3ScheduleKey(item.schedule)),state.fresh_queue))
   state.fresh_queue=[{item} state.fresh_queue];
  end
 end
 state.construction_history{end}.queued_candidates=construction.queued_candidates;
 li=find(strcmp(state.lineage_keys,node.schedule.root_key),1);
 if isempty(li)
  state.lineage_keys{end+1}=node.schedule.root_key; li=numel(state.lineage_keys);
  state.lineage_visits(li)=-1; state.lineage_stalls(li)=0;
 end
 if construction.consistency.passed&&sum(node.visited)>state.lineage_visits(li)
  state.lineage_visits(li)=sum(node.visited); state.lineage_stalls(li)=0;
 else
  state.lineage_stalls(li)=state.lineage_stalls(li)+1;
 end
 if state.lineage_stalls(li)>=c.lineage_stall_limit
  if c.replan_enabled&&sum(node.visited)>=2&&~all(node.visited)
   enqueue(node.schedule,1:35,[],c.recovery_stall_limit);
  end
  keep=~cellfun(@(n)strcmp(n.schedule.root_key,node.schedule.root_key),beam); beam=beam(keep);
  if isempty(beam)
   beam=ctocscreen.v3Roots(eph,c,stream); state.root_archive=[state.root_archive beam];
   state.pending_roots=[state.pending_roots beam];
  end
 end
 if construction.consistency.passed && (isempty(state.best_partial)||sum(node.visited)>sum(state.best_partial.visited) ...
    ||(sum(node.visited)==sum(state.best_partial.visited)&&node.J<state.best_partial.J))
  state.best_partial=node;
 end
 for feedbackIndex=1:numel(constructed)
  [pheromone,state.prefix_population,prefixFeedback]=ctocscreen.v3PrefixReinforce( ...
   pheromone,state.prefix_population,constructed{feedbackIndex},eph,c);
  prefixFeedback.iteration=iter; state.prefix_feedback_history{end+1}=prefixFeedback;
 end
 % A complete fresh candidate is ready only when no rejected recovery task is
 % waiting. Otherwise fresh FIFO must not starve the full-task replan queue.
 completeReady=isempty(state.work_pool)&&any(cellfun(@(p)numel(p.target_ids)==35,state.fresh_queue)) ...
  &&isempty(state.recovery_pool)&&isempty(state.replan_queue);
 [state,source]=ctocscreen.v3BeginAttempt(state,~isempty(state.work_pool)&&~completeReady, ...
  (~isempty(state.recovery_pool)||(c.replan_enabled&&~isempty(state.replan_queue)))&&~completeReady, ...
  ~isempty(state.pending_initial)&&~completeReady);
 attemptId=numel(state.source_history); recoveryStalls=0; taskIds=1:35; queuedFresh=false;
 try
 if strcmp(source,'initial')
  current=[]; seed=state.pending_initial{1}; state.pending_initial(1)=[];
  operator='initial_candidate'; proposalKeys={};
 elseif strcmp(source,'work')
  current=state.work_pool{randi(stream,numel(state.work_pool))};
  [state,workOperator]=ctocscreen.v3WorkOperator(state,current.schedule,c);
  if strcmp(workOperator,'joint_refine')
   seed=current.schedule; operator='joint_refine';
  elseif strcmp(workOperator,'directed_replan')
   [rebuilt,~]=rebuild(current.schedule,'work',current,1:35);
   if isempty(rebuilt)
    outcome('replan_failed'); colonies{active}=beam; checkpoint(); continue;
   end
   seed=rebuilt{1}.schedule; operator=['directed_' rebuilt{1}.mode];
   for j=2:numel(rebuilt), enqueue(rebuilt{j}.schedule,1:35,current,0); end
  else
   [seed,operator,connection]=ctocscreen.v3ConnectedMutation(current.schedule,eph,c,stream,min(c.joint_seconds,c.budget_s-toc(clock)));
   state.mutation_history{end+1}=struct('operator',operator,'connection',connection);
   if ~connection.passed
    if connection.difference.novel||connection.visit_assignment_changed, enqueue(seed,taskIds,current,0); end
    outcome('construction_failed'); colonies{active}=beam; checkpoint(); continue;
   end
  end
  proposalKeys={};
 elseif strcmp(source,'recovery')
  recoveryAttempt=sum(cellfun(@(h)strcmp(h.source,'recovery'),state.source_history));
  if c.replan_enabled&&~isempty(state.replan_queue)&&(isempty(state.recovery_pool)||mod(recoveryAttempt,2)==0)
   recovery=state.replan_queue{1}; state.replan_queue(1)=[]; current=recovery.parent;
   [rebuilt,replan]=rebuild(recovery.schedule,'stalled_recovery',current,recovery.target_ids);
   if isempty(rebuilt)
    if replan.attempt_number<c.replan_retry_limit
     [state.replan_queue,~]=ctocscreen.v3QueueCandidates(state.replan_queue,{recovery},c.archive_size);
    else
     state.replanned_keys{end+1}=replan.task_key;
    end
    outcome('replan_failed'); colonies{active}=beam; checkpoint(); continue;
   end
   state.replanned_keys{end+1}=replan.task_key;
   seed=rebuilt{1}.schedule; taskIds=1:35; operator=['directed_' rebuilt{1}.mode];
   for j=1:numel(rebuilt)
    if isfield(rebuilt{j},'node')
     % Proposal nodes are not certified prefixes. They may seed construction
     % and must pass the usual consistency/full-task gates before credit.
     beam=ctocscreen.v3SelectBeam([beam {rebuilt{j}.node}],c.beam_width,c);
    end
   end
   for j=2:numel(rebuilt), enqueue(rebuilt{j}.schedule,1:35,current,0); end
  else
   slot=randi(stream,numel(state.recovery_pool));
   recovery=state.recovery_pool{slot}; state.recovery_pool(slot)=[];
   seed=recovery.schedule; current=recovery.parent; taskIds=recovery.target_ids;
   recoveryStalls=recovery.stalls; operator='resume_recovery';
  end
  proposalKeys={};
 else
  current=[];
  if isempty(state.fresh_queue)
   seed=node.schedule; operator='beam_completion'; taskIds=find(isfinite(seed.witness_times_s)).';
  else
   slot=find(cellfun(@(p)numel(p.target_ids)==35,state.fresh_queue),1);
   if isempty(slot), slot=1; end
   fresh=state.fresh_queue{slot}; state.fresh_queue(slot)=[];
   queuedFresh=true;
   seed=fresh.schedule; node=fresh.node; operator=fresh.operator; taskIds=fresh.target_ids;
  end
  proposalKeys=node.keys;
 end
 catch err
  state.failures{end+1}=struct('id',err.identifier,'message',err.message,'stage','mutation');
  outcome('proposal_failed'); colonies{active}=beam; checkpoint(); continue;
 end
 if toc(clock)>=c.budget_s
  pending=struct('schedule',seed,'target_ids',taskIds,'parent',current,'stalls',recoveryStalls);
  if strcmp(source,'fresh')
   pending=struct('schedule',seed,'target_ids',taskIds,'node',node,'operator',operator);
  end
  state=ctocscreen.v3DeferProposal(state,source,pending);
  outcome('budget_exhausted_deferred'); break;
 end
 cc=c; cc.joint_seconds=min(c.joint_seconds,c.budget_s-toc(clock));
 cc.joint_target_ids=taskIds;
 if sum(vecnorm(seed.delta_v_km_s,2,2))>c.search_max_dv_km_s&& ...
    numel(taskIds)<35&&~strcmp(source,'recovery')
  isNew=rejectSchedule(seed,source);
  if isNew, enqueue(seed,1:35,current,c.recovery_stall_limit); end
  outcome('user_cost_rejected_to_recovery'); colonies{active}=beam; checkpoint(); continue;
 end
 if c.cold_completion_enabled&&numel(taskIds)==35&&sum(isfinite(seed.witness_times_s))<35
  % Restore the missing trajectory structure before solving a 35-target NLP.
  % A few visited targets are not a physically initialized full mission.
  try
   tail=ctocscreen.v3PrefixNode(seed,seed.duration_s,eph,c);
   tailConfig=c; tailConfig.completion_seconds=min(c.completion_seconds,max(.001,c.budget_s-toc(clock)));
   [tail,tailReport,tailFailures,tailNodes]=ctocscreen.v3Complete(tail,eph,tailConfig,stream,pheromone,true,true);
   tailReport.iteration=iter; tailReport.target_ids=1:35;
   state.full_construction_history{end+1}=tailReport; state.failures=[state.failures tailFailures];
   if tailReport.consistency.passed
    seed=tail.schedule; beam=ctocscreen.v3SelectBeam([tailNodes beam],c.beam_width,c);
    if isempty(state.best_partial)||sum(tail.visited)>sum(state.best_partial.visited), state.best_partial=tail; end
   end
   if ~tailReport.consistency.passed||~all(tail.visited)
    enqueue(seed,1:35,current,recoveryStalls+1);
    outcome('full_construction_pending'); colonies{active}=beam; checkpoint(); continue;
   end
  catch err
   state.failures{end+1}=struct('id',err.identifier,'message',err.message,'stage','full_construction');
   enqueue(seed,1:35,current,c.recovery_stall_limit);
   outcome('full_construction_failed'); colonies{active}=beam; checkpoint(); continue;
  end
  cc.joint_seconds=min(c.joint_seconds,max(0,c.budget_s-toc(clock)));
 end
 % Partial prefixes only receive construction refinement. The complete NLP
 % always includes all 35 targets, and SA/export still require true 35/35.
 if numel(taskIds)<35
  if numel(taskIds)<2 || (strcmp(source,'fresh')&& ...
     mod(state.source_history{attemptId}.source_attempt,c.prefix_joint_every)~=0)
   if queuedFresh&&numel(taskIds)>=2, state.fresh_queue{end+1}=fresh; end
   outcome('prefix_deferred');
   colonies{active}=beam; checkpoint(); continue
  end
  cc.joint_seconds=min(cc.joint_seconds,c.prefix_joint_seconds);
  operator=['prefix_' operator];
 end
 inputKey=ctocscreen.v3ScheduleKey(seed);
 if toc(clock)>=c.budget_s
  if strcmp(source,'fresh')
   pending=struct('schedule',seed,'target_ids',taskIds,'node',node,'operator',operator);
  else
   pending=struct('schedule',seed,'target_ids',taskIds,'parent',current,'stalls',recoveryStalls);
  end
  state=ctocscreen.v3DeferProposal(state,source,pending);
  outcome('budget_exhausted_deferred'); break;
 end
 if strcmp(source,'fresh')&&any(strcmp(state.refined_keys,inputKey))
  state.duplicate_joint_skips=state.duplicate_joint_skips+1;
  outcome('duplicate_input');
  colonies{active}=beam; checkpoint(); continue
 end
 constructedCheck=[];
 if isempty(state.work_pool)&&numel(taskIds)==35&&all(isfinite(seed.witness_times_s))&& ...
    ~isfield(seed,'visit_plan_times_s')&&sum(vecnorm(seed.delta_v_km_s,2,2))<=c.search_max_dv_km_s
  verifyClock=tic; constructedCheck=ctocscreen.v3Verify(seed,eph,c);
 end
 if ~isempty(constructedCheck)&&constructedCheck.passed
  cand=seed;
  diagnostic=struct('status','screened_feasible','refinement_performed',false, ...
   'selected_source','constructed_complete','elapsed_s',toc(verifyClock),'replay',constructedCheck, ...
   'active_target_ids',taskIds,'verification_only',true);
 else
  cc.joint_seconds=min(cc.joint_seconds,max(0,c.budget_s-toc(clock)));
  [cand,diagnostic]=ctocscreen.v3JointOptimize(seed,eph,cc);
 end
 if diagnostic.refinement_performed&&~any(strcmp(state.refined_keys,inputKey))
  state.refined_keys{end+1}=inputKey;
 end
 record=struct('iteration',iter,'colony',active,'operator',operator,'diagnostic',diagnostic, ...
  'accepted',false,'sa_probability',NaN,'sa_draw',NaN,'temperature',state.temperature, ...
  'exported',false,'verification',[],'proposal_keys',{proposalKeys});
 record.source_family=source; record.source_attempt_id=attemptId;
 record.input_schedule_key=inputKey; record.task_target_ids=taskIds;
 record.from_candidate_queue=queuedFresh;
 if isfield(seed,'replan_origin'), record.replan_origin=seed.replan_origin; end
 record.parent_J_km_s=NaN;
 if ~isempty(current), record.parent_J_km_s=current.J; end
 if isfield(diagnostic,'verification_only')&&diagnostic.verification_only, outcome('construction_verified');
 else, outcome('joint_completed'); end
 if sum(vecnorm(cand.delta_v_km_s,2,2))>c.search_max_dv_km_s
  rejectSchedule(cand,'joint_output'); record.policy_rejected=true;
  if strcmp(diagnostic.status,'screened_feasible')
   checked=ctocscreen.v3Verify(cand,eph,c); record.verification=checked;
   if checked.passed&&(isempty(state.best_rejected_complete)||checked.total_dv_km_s<state.best_rejected_complete.verification.total_dv_km_s)
    state.best_rejected_complete=struct('schedule',cand,'verification',checked, ...
     'elapsed_s',baseElapsed+toc(clock),'accepted',false,'reason','user_cost_threshold');
    rejected=state.best_rejected_complete; save(fullfile(folder,'rejected_complete.mat'),'rejected','-v7.3');
   end
  end
  enqueue(cand,1:35,current,c.recovery_stall_limit);
  outcome('user_cost_rejected'); state.history{end+1}=record; colonies{active}=beam; checkpoint(); continue;
 end
 if strcmp(diagnostic.status,'screened_feasible')
  % A screening-only complete schedule must not seed SA or mutations before
  % fixed-pulse independent verification. Rejected schedules stay diagnostic.
  if ~isempty(constructedCheck)&&constructedCheck.passed, verified=constructedCheck;
  else, verified=ctocscreen.v3Verify(cand,eph,c); end
  record.verification=verified;
  if ~verified.passed
   enqueue(cand,taskIds,current,recoveryStalls+1); outcome('independent_rejected');
   record.independent_rejected=true; state.history{end+1}=record;
   colonies{active}=beam; checkpoint(); continue
  end
  if ~isempty(current)
   record.physical_difference=ctocscreen.v3MutationDifference(current.schedule,cand,c,eph.model);
   if ~record.physical_difference.novel && ~strcmp(operator,'joint_refine')
    record.retained_existing=true; state.history{end+1}=record; outcome('returned_to_parent');
    colonies{active}=beam; checkpoint(); continue
   end
  end
  if isfield(cand,'visit_plan_times_s'), cand=rmfield(cand,'visit_plan_times_s'); end
  cand.witness_times_s=verified.witness_times_s;
  J=verified.total_dv_km_s;
  if isempty(current), accepted=true;
  else
   [accepted,record.sa_probability,record.sa_draw]=ctocscreen.v3SAAccept(current.J,J,state.temperature,stream);
   if abs(J-current.J)>1e-10, state.temperature_samples(end+1)=abs(J-current.J); end
   if numel(state.temperature_samples)<=5 && any(state.temperature_samples>0)
    state.temperature=max(1e-6,median(state.temperature_samples(state.temperature_samples>0)));
   end
  end
  candidateKey=ctocscreen.v3ScheduleKey(cand);
  novel=~any(cellfun(@(w)strcmp(ctocscreen.v3ScheduleKey(w.schedule),candidateKey),state.work_pool));
  record.retained_existing=~novel;
  record.sa_accepted=accepted;
  item=deposit(cand,J); record.feedback=state.feedback_history{end};
  if accepted&&novel
   [state.work_pool,record.pool_admission]=ctocscreen.v3SelectWorkPool(state.work_pool,item,c,~isempty(current));
   state.pool_history{end+1}=record.pool_admission;
   record.accepted=record.pool_admission.retained;
  end
  if isempty(state.elite)||J<state.elite.verification.total_dv_km_s
   if verified.passed
    state.elite=struct('schedule',cand,'verification',verified,'elapsed_s',baseElapsed+toc(clock));
    state.export_archive{end+1}=state.elite;
    elite=state.elite; save(fullfile(folder,'elite.pending.mat'),'elite','-v7.3');
    movefile(fullfile(folder,'elite.pending.mat'),fullfile(folder,'elite.mat'),'f');
    record.exported=true; state.stalled=0;
   end
  end
 else
  if isfield(diagnostic,'recovery_progress')&&diagnostic.recovery_progress, recoveryStalls=0;
  else, recoveryStalls=recoveryStalls+1; end
  if numel(taskIds)==35||~strcmp(diagnostic.status,'prefix_restored')
   enqueue(cand,taskIds,current,recoveryStalls);
  end
  if numel(taskIds)<35 && ~isempty(diagnostic.replay)&&diagnostic.replay.height_passed&&strcmp(diagnostic.selected_source,'solver_candidate')
   repaired=node; repaired.schedule=cand; repaired.state=diagnostic.replay.final_state;
   repaired.t=cand.duration_s; repaired.visited=diagnostic.replay.distance_km<=1;
   repaired.J=diagnostic.replay.total_dv_km_s; repaired.keys={};
   repaired.inclination_penalty=diagnostic.replay.inclination_penalty;
   repaired.estimate=repaired.J+c.plane_penalty_km_s*repaired.inclination_penalty; repaired.last_gain=0;
   if strcmp(source,'fresh')&&strcmp(diagnostic.selected_source,'solver_candidate')
    beam=ctocscreen.v3ReplacePrefix(beam,seed,repaired);
   else
    beam=[beam,{repaired}];
   end
   beam=ctocscreen.v3SelectBeam(beam,c.beam_width,c);
  end
 end
 state.history{end+1}=record; state.stalled=state.stalled+1;
 state.temperature=max(1e-6,state.temperature*c.cooling);
 if mod(state.stalled,c.reheat_after)==0, state.temperature=max(c.temperature_km_s,2*state.temperature); end
 % Migration: new geometry returns to construction through roots of work pool.
 if ~isempty(state.work_pool)&&mod(iter,3)==0
  root=ctocscreen.v3Roots(eph,ctocscreen.v3Defaults(struct('root_count',1)),stream); root=root{1};
  migrant=randi(stream,numel(state.work_pool));
  root.schedule.initial_q=state.work_pool{migrant}.schedule.initial_q;
  root.schedule.root_key=sprintf('%.17g/',root.schedule.initial_q);
  root.seed_target=0; root.seed_duration_s=NaN; root.root_method='work_orbit_reuse';
  root.state=ctocscreen.initialState(root.schedule.initial_q,eph.model.mu,eph.model.re);
  beam=ctocscreen.v3SelectBeam([beam {root}],c.beam_width,c);
 end
 state.config=c;
 colonies{active}=beam;
 if mod(iter,c.checkpoint_every)==0, checkpoint(); end
end
checkpoint();
fid=fopen(fullfile(folder,'summary.md'),'w','n','UTF-8'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# V3 search run\n\nDataset: %s. Nominal J2; official ATK alignment pending.\n\n',state.dataset_kind);
fprintf(fid,'Seed %d; iterations %d; cumulative wall %.3f s; current budget overrun %.3f s.\n\n', ...
 c.seed,state.iteration,state.elapsed_s,state.budget_overrun_s);
jointCalls=sum(cellfun(@(h)~isfield(h.diagnostic,'verification_only')||~h.diagnostic.verification_only,state.history));
fprintf(fid,'Candidate records: %d; joint calls: %d; retained work %d, recovery %d; failures recorded %d.\n\n', ...
 numel(state.history),jointCalls,numel(state.work_pool),numel(state.recovery_pool),numel(state.failures));
fprintf(fid,'Source attempts: %d (failed/skipped attempts advance scheduling).\n\n',numel(state.source_history));
fprintf(fid,'Pending constructed candidates: %d; directed replans: %d; pending replan requests: %d.\n\n', ...
 numel(state.fresh_queue),numel(state.replan_history),numel(state.replan_queue));
if isempty(state.elite)
 fprintf(fid,'No independently verified 35-target elite. No elite.mat published.\n');
else
 fprintf(fid,'Independent nominal J2 elite: 35/35, %.12g km/s. Not ATK verified.\n',state.elite.verification.total_dv_km_s);
end
 function item=deposit(schedule,J)
  [experience,features,actions]=ctocscreen.v3ScheduleExperience(schedule,eph,c);
  [pheromone,state.reinforcement_ledger,feedback]=ctocscreen.v3Reinforce( ...
   pheromone,state.reinforcement_ledger,features,experience,J,c);
  feedback.iteration=state.iteration;
  feedback.plane_policy={};
  for j=1:numel(actions)
   a=actions{j};
   feedback.plane_policy{end+1}=ctocscreen.v3PolicyFeedback(pheromone,a.origin,a.target,a.dt,a.dv,c);
  end
  state.feedback_history{end+1}=feedback;
  item=struct('schedule',schedule,'J',J,'features',features);
 end
 function enqueue(schedule,ids,parent,stalls)
  state=ctocscreen.v3QueueRecovery(state,schedule,ids,parent,stalls,c);
 end
 function isNew=rejectSchedule(schedule,stage,checkedPrefix)
  raw=sum(vecnorm(schedule.delta_v_km_s,2,2));
  key=ctocscreen.v3ScheduleKey(schedule); isNew=~any(strcmp(state.policy_rejected_keys,key));
  if ~isNew, return; end
  state.policy_rejected_keys{end+1}=key;
  item=struct('schedule',schedule,'total_dv_km_s',raw,'stage',stage, ...
   'reason','user_search_cost_threshold','physical_feasibility','not_certified_by_rejection', ...
   'feedback',{{}});
  try
   reusable=nargin>=3&&checkedPrefix.passed&&isfield(checkedPrefix,'policy_actions')&& ...
    strcmp(checkedPrefix.schedule_key,key)&&strcmp(checkedPrefix.witness_key,mat2str(schedule.witness_times_s,17));
   item.reused_prefix_states=reusable;
   if reusable
    actions=checkedPrefix.policy_actions; complete=all(isfinite(schedule.witness_times_s));
   else
    [replay,trace]=ctocscreen.v3Replay(schedule,eph,c,false,'prefix_witnesses');
    assert(~strcmp(replay.status,'propagation_failure'),'ctocscreen:v3:feedbackReplay','No negative evidence from failed propagation.');
    [~,actions]=ctocscreen.v3TraceExperience(schedule,trace,c); complete=replay.visit_count==35;
   end
   for j=1:numel(actions)
    a=actions{j};
    item.feedback{end+1}=ctocscreen.v3PolicyFeedback(pheromone,a.origin,a.target,a.dt,a.dv,c,raw,complete);
   end
  catch err
   item.feedback_failure=err.identifier;
  end
  state.policy_rejection_count=state.policy_rejection_count+isNew;
  state.policy_rejections{end+1}=item;
  if numel(state.policy_rejections)>c.archive_size, state.policy_rejections(1)=[]; end
 end
 function [proposals,report]=rebuild(schedule,family,parent,ids)
  key=ctocscreen.v3RecoveryKey(schedule,ids,parent);
  attempt=1+sum(cellfun(@(h)strcmp(h.task_key,key),state.replan_history));
  cutRound=1+sum(cellfun(@(h)isfield(h,'root_key')&&strcmp(h.root_key,schedule.root_key),state.replan_history));
  [proposals,report]=ctocscreen.v3Replan(schedule,eph,c,stream,pheromone, ...
   min(c.replan_seconds,max(0,c.budget_s-toc(clock))),attempt,cutRound);
  report.task_key=key;
  report.iteration=iter; report.source=family; report.parent_J_km_s=NaN;
  report.parent_key='';
  if ~isempty(parent)
   report.parent_J_km_s=parent.J; report.parent_key=ctocscreen.v3ScheduleKey(parent.schedule);
  end
  state.replan_history{end+1}=report;
 end
 function outcome(value)
  state.source_history{attemptId}.outcome=value;
 end
 function checkpoint()
  state.config=c;
  state.rng_state=stream.State; state.beam=beam;
  colonies{active}=beam; state.colonies=colonies;
  state.pheromone_keys=keys(pheromone); state.pheromone_values=values(pheromone);
  state.elapsed_s=baseElapsed+toc(clock); state.budget_overrun_s=max(0,toc(clock)-c.budget_s);
  save(fullfile(folder,'checkpoint.pending.mat'),'state','-v7.3');
  movefile(fullfile(folder,'checkpoint.pending.mat'),fullfile(folder,'checkpoint.mat'),'f');
 end
end
