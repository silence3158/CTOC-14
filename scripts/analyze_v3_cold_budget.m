function report=analyze_v3_cold_budget(label,comparisonLabel)
%ANALYZE_V3_COLD_BUDGET Diagnose saved runs; never seed or execute a search.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
folder=fullfile(sim,'runs/v3/search',label); data=load(fullfile(folder,'checkpoint.mat'),'state'); s=data.state;
report=struct('purpose','saved_run_diagnosis_only','label',label,'elapsed_s',s.elapsed_s);
report.first_complete_s=s.export_archive{1}.elapsed_s;
report.first_J=s.export_archive{1}.verification.total_dv_km_s;
report.best_J=s.elite.verification.total_dv_km_s;
report.improvement_km_s=report.first_J-report.best_J;
report.improvement_percent=100*report.improvement_km_s/report.first_J;
hist=s.history; work=cellfun(@(h)strcmp(h.source_family,'work'),hist);
refine=cellfun(@(h)strcmp(h.operator,'joint_refine'),hist);
report.work_calls=sum(work); report.work_refinements=sum(refine);
report.refinement_s=sum(cellfun(@(h)h.diagnostic.elapsed_s,hist(refine)));
report.joint_record_s=sum(cellfun(@(h)h.diagnostic.elapsed_s,hist));
report.sa_trials=sum(cellfun(@(h)isfinite(h.sa_probability),hist));
report.sa_uphill_trials=sum(cellfun(@(h)isfinite(h.sa_probability)&& ...
 ~isempty(h.verification)&&h.verification.total_dv_km_s>h.parent_J_km_s,hist));
report.mutations=numel(s.mutation_history); report.replans=numel(s.replan_history);
report.new_complete_replans=sum(cellfun(@(h)startsWith(h.operator,'directed_')&&h.accepted,hist));
report.root_methods=cellfun(@(h)h.root_method,s.construction_history,'UniformOutput',false);
report.construction=zeros(numel(s.construction_history),7);
for k=1:numel(s.construction_history)
 h=s.construction_history{k};
 report.construction(k,:)=[h.iteration h.protected_root h.coverage_completion h.start_visits h.end_visits h.elapsed_s h.started_at_s];
 fprintf('CONSTRUCT it%2d protected%d coverageMode%d %2d->%2d wall%.3f root%s\n', ...
  report.construction(k,1:6),h.root_method);
end
orders=zeros(numel(s.work_pool),35); roots=cell(1,numel(s.work_pool));
report.pool=zeros(numel(s.work_pool),5);
for k=1:numel(s.work_pool)
 w=s.work_pool{k}; [~,order]=sort(w.schedule.witness_times_s); orders(k,:)=order.';
 roots{k}=w.schedule.root_key;
 report.pool(k,:)=[w.J,numel(w.schedule.maneuver_times_s),w.schedule.duration_s, ...
  max(vecnorm(w.schedule.delta_v_km_s,2,2)),norm(w.schedule.initial_q-s.export_archive{1}.schedule.initial_q)];
 fprintf('POOL %d J%.9f burns%d T%.1f maxDv%.4f root%s\n',k,report.pool(k,1:4),roots{k});
end
report.unique_pool_orders=size(unique(orders,'rows'),1); report.unique_pool_lineages=numel(unique(roots));
report.reward_reasons=cellfun(@(h)h.reason,s.feedback_history,'UniformOutput',false);
report.mutation_details=s.mutation_history;
if nargin>=2
 b=load(fullfile(sim,'runs/v3/search',comparisonLabel,'checkpoint.mat'),'state'); b=b.state;
 report.comparison=struct('label',comparisonLabel,'same_signature',isequaln(s.signature,b.signature), ...
  'same_initial_roots',isequaln(s.root_archive(1:8),b.root_archive(1:8)),'config_differences',{{}});
 names=fieldnames(s.config);
 for k=1:numel(names)
  name=names{k};
  if ~isequaln(s.config.(name),b.config.(name))
   report.comparison.config_differences{end+1}=struct('field',name,'run',s.config.(name),'comparison',b.config.(name));
   fprintf('CONFIG_DIFF %s\n',name);
  end
 end
 report.comparison.first_pairs_equal=true; report.comparison.first_pair_difference=[];
 a=s.construction_history{1}.expansion_history; bb=b.construction_history{1}.expansion_history;
 for k=1:min(numel(a),numel(bb))
  p=a{k}.time_search{1}.selected_pairs; q=bb{k}.time_search{1}.selected_pairs;
  pa=cellfun(@(v)[v.target v.dt],p,'UniformOutput',false);
  pb=cellfun(@(v)[v.target v.dt],q,'UniformOutput',false);
  if ~isequaln(pa,pb)
   report.comparison.first_pairs_equal=false;
   report.comparison.first_pair_difference=struct('expansion',k,'run',{pa},'comparison',{pb});
   fprintf('FIRST_PAIR_DIFFERENCE expansion%d\n',k); break;
  end
 end
 fprintf('COMPARE sameSource%d sameRoots%d firstPairsEqual%d\n', ...
  report.comparison.same_signature,report.comparison.same_initial_roots,report.comparison.first_pairs_equal);
end
fprintf('SUMMARY first%.12f best%.12f gain%.9f percent%.6f work%d refine%d refineWall%.3f totalRecordsWall%.3f\n', ...
 report.first_J,report.best_J,report.improvement_km_s,report.improvement_percent, ...
 report.work_calls,report.work_refinements,report.refinement_s,report.joint_record_s);
fprintf('DIVERSITY orders%d lineages%d mutations%d replans%d acceptedReplan%d SA%d uphillSA%d\n', ...
 report.unique_pool_orders,report.unique_pool_lineages,report.mutations,report.replans, ...
 report.new_complete_replans,report.sa_trials,report.sa_uphill_trials);
save(fullfile(folder,'budget_diagnosis.mat'),'report');
end
