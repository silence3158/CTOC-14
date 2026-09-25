function report=report_v3_cold_run(label)
%REPORT_V3_COLD_RUN Report full-mission success separately from construction.
sim=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(sim,'runs','v3','search',label); data=load(fullfile(folder,'checkpoint.mat'),'state'); s=data.state;
report=struct('seed',s.config.seed,'cold_start',isempty(s.config.initial_candidates)&&isempty(s.config.resume_file), ...
 'elapsed_s',s.elapsed_s,'iterations',s.iteration,'complete_passed',false,'first_complete_s',NaN, ...
 'first_complete_dv_km_s',NaN, ...
 'total_dv_km_s',NaN,'best_prefix_visits',0,'prefix_joint_calls',0,'full_joint_calls',0,'verification_only',0, ...
 'construction_s',0,'prefix_check_s',0,'prefix_checks',0,'prefix_cache_hits',0,'truncated_checks',0, ...
 'prefix_admissions',0,'full_feedback_events',numel(s.feedback_history), ...
 'full_rewards',sum(cellfun(@(h)h.amount_per_key>0,s.feedback_history)),'pheromone_queries',0, ...
 'nonneutral_queries',0,'prefix_guidance_hits',0,'mutations',numel(s.mutation_history));
report.negative_guidance_hits=0; report.actual_cost_rejections=0; report.actual_plane_penalties=0;
report.new_negative_feedback=0;
report.rejected_complete_verified=false; report.rejected_complete_dv_km_s=NaN;
if isfield(s,'best_rejected_complete')&&~isempty(s.best_rejected_complete)
 report.rejected_complete_verified=s.best_rejected_complete.verification.passed;
 report.rejected_complete_dv_km_s=s.best_rejected_complete.verification.total_dv_km_s;
end
report.search_max_dv_km_s=Inf;
if isfield(s.config,'search_max_dv_km_s'), report.search_max_dv_km_s=s.config.search_max_dv_km_s; end
if ~isempty(s.elite)
 report.complete_passed=s.elite.verification.passed;
 report.total_dv_km_s=s.elite.verification.total_dv_km_s;
 report.first_complete_s=s.export_archive{1}.elapsed_s;
 report.first_complete_dv_km_s=s.export_archive{1}.verification.total_dv_km_s;
end
if ~isempty(s.best_partial), report.best_prefix_visits=sum(s.best_partial.visited); end
onlyVerify=cellfun(@(h)isfield(h.diagnostic,'verification_only')&&h.diagnostic.verification_only,s.history);
report.verification_only=sum(onlyVerify);
report.prefix_joint_calls=sum(cellfun(@(h)numel(h.task_target_ids)<35,s.history(~onlyVerify)));
report.full_joint_calls=sum(~onlyVerify)-report.prefix_joint_calls;
if isfield(s,'prefix_feedback_history')
 report.prefix_admissions=sum(cellfun(@(h)h.admitted,s.prefix_feedback_history));
end
histories=s.construction_history;
if isfield(s,'full_construction_history'), histories=[histories s.full_construction_history]; end
for k=1:numel(s.replan_history)
 for j=1:numel(s.replan_history{k}.attempts)
  completion=s.replan_history{k}.attempts{j}.completion;
  if ~isempty(completion), histories{end+1}=completion; end
 end
end
report.prefix_restarts=0;
if isfield(s,'prefix_restart_history'), report.prefix_restarts=numel(s.prefix_restart_history); end
for k=1:numel(histories)
 h=histories{k}; report.construction_s=report.construction_s+h.elapsed_s;
 report.prefix_checks=report.prefix_checks+h.consistency_checks;
 if isfield(h,'prefix_check_seconds'), report.prefix_check_s=report.prefix_check_s+h.prefix_check_seconds; end
 if isfield(h,'consistency_cache_hits'), report.prefix_cache_hits=report.prefix_cache_hits+h.consistency_cache_hits; end
 report.truncated_checks=report.truncated_checks+sum(cellfun(@(r)r.truncated,h.candidate_checks));
 for j=1:numel(h.expansion_history), addQueries(h.expansion_history{j}); end
end
for k=1:numel(s.beam_expansion_history), addQueries(s.beam_expansion_history{k}); end
if isfield(s,'policy_rejections')
 for k=1:numel(s.policy_rejections), addFeedback(s.policy_rejections{k}.feedback); end
end
for k=1:numel(s.feedback_history)
 if isfield(s.feedback_history{k},'plane_policy'), addFeedback(s.feedback_history{k}.plane_policy); end
end
fprintf('Cold=%d seed=%d wall=%.3f s iterations=%d complete=%d prefix=%d/35\n', ...
 report.cold_start,report.seed,report.elapsed_s,report.iterations,report.complete_passed,report.best_prefix_visits);
fprintf('First complete=%.3f s DV=%.9f km/s; best DV=%.9f km/s; joint full=%d prefix=%d\n', ...
 report.first_complete_s,report.first_complete_dv_km_s,report.total_dv_km_s,report.full_joint_calls,report.prefix_joint_calls);
fprintf('Construction=%.3f s checks=%.3f s n=%d cacheHits=%d truncated=%d\n', ...
 report.construction_s,report.prefix_check_s,report.prefix_checks,report.prefix_cache_hits,report.truncated_checks);
fprintf('Prefix admissions=%d full rewards=%d queries=%d nonneutral=%d prefixHits=%d mutations=%d\n', ...
 report.prefix_admissions,report.full_rewards,report.pheromone_queries,report.nonneutral_queries,report.prefix_guidance_hits,report.mutations);
fprintf('Policy cap=%.6g km/s cost-penalized action observations=%d plane penalties=%d negative updates=%d negative hits=%d\n', ...
 report.search_max_dv_km_s,report.actual_cost_rejections,report.actual_plane_penalties,report.new_negative_feedback,report.negative_guidance_hits);
if ~isempty(s.best_partial), fprintf('Construction archive (not acceptance): DV %.9g km/s duration %.3f s\n',s.best_partial.J,s.best_partial.t); end
fprintf('Independently verified but policy-rejected complete: %d, DV %.9g km/s (not an elite).\n', ...
 report.rejected_complete_verified,report.rejected_complete_dv_km_s);
save(fullfile(folder,'cold_report.mat'),'report');

 function addQueries(expansion)
  if isfield(expansion,'policy_feedback'), addFeedback(expansion.policy_feedback); end
  for a=1:numel(expansion.time_search)
   screen=expansion.time_search{a}.screen;
   if isfield(screen,'pheromone_queries')
    report.pheromone_queries=report.pheromone_queries+screen.pheromone_queries;
    report.nonneutral_queries=report.nonneutral_queries+screen.nonneutral_queries;
    report.prefix_guidance_hits=report.prefix_guidance_hits+screen.prefix_guidance_hits;
    if isfield(screen,'negative_guidance_hits'), report.negative_guidance_hits=report.negative_guidance_hits+screen.negative_guidance_hits; end
   end
   search=expansion.time_search{a};
   if isfield(search,'policy_feedback')
    addFeedback(search.policy_feedback);
   end
  end
 end
 function addFeedback(feedback)
  for j=1:numel(feedback)
   f=feedback{j};
   report.actual_cost_rejections=report.actual_cost_rejections+f.rejected_cost;
   report.actual_plane_penalties=report.actual_plane_penalties+(f.inclination_change_deg>s.config.plane_change_threshold_deg);
   report.new_negative_feedback=report.new_negative_feedback+(f.new_penalty>0);
  end
 end
end
