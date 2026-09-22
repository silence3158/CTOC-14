function report=report_structure_experiment(folder)
%REPORT_STRUCTURE_EXPERIMENT Pair replacements against unchanged-topology R.
d=load(fullfile(folder,'design.mat'));W=height(d.manifest);arms=d.cfg.arms;
rows=cell(W*numel(arms),18);cost=nan(W,numel(arms));available=false(W,numel(arms));at=0;
for b=1:numel(d.sources)
 writetable(ctocscreen.burnDiagnostics(d.sources{b},d.problems{b}),fullfile(folder,sprintf('source%02d_burns.csv',b)));
end
for w=1:W
 b=d.manifest.source(w);baseline=d.sources{b}.independent.total_dv_km_s;
 for a=1:numel(arms)
  at=at+1;status='not_completed';dv=NaN;budget=NaN;rep=NaN;feasible=0;elapsed=NaN;parts=[NaN NaN NaN];checks=nan(1,5);
  file=fullfile(folder,sprintf('window%02d_%s',w,arms(a)),'result.mat');
  if isfile(file)
   x=load(file,'out');o=x.out;status=o.status;dv=o.best.independent.total_dv_km_s;
   budget=o.budget_best.independent.total_dv_km_s;feasible=o.feasible_replacements;elapsed=o.elapsed_s;
   if ~isempty(o.budget_replacement)
    br=o.budget_replacement;rep=br.independent.total_dv_km_s;parts=br.cost_parts_km_s;available(w,a)=true;
    checks=[br.independent.unique_visit_count max(br.independent.event_distances_km)*1000 br.independent.min_altitude_km br.schedule.duration_s/86400 numel(br.schedule.maneuver_times_s)];
    writetable(ctocscreen.burnDiagnostics(br,d.problems{b}),fullfile(fileparts(file),'burns_budget_replacement.csv'));
   end
   if ~isempty(o.best_replacement)
    writetable(ctocscreen.burnDiagnostics(o.best_replacement,d.problems{b}),fullfile(fileparts(file),'burns_final_replacement.csv'));
   end
   if strcmp(status,'completed'),cost(w,a)=budget;end
  end
  rows(at,:)=[{w,b,arms(a),status,baseline,budget,dv,rep,feasible,elapsed,parts(1),parts(2),parts(3)} num2cell(checks)];
 end
end
tasks=cell2table(rows,'VariableNames',{'window','source','arm','status','baseline_km_s','budget_best_km_s','final_best_km_s', ...
 'budget_replacement_km_s','feasible_replacement_trials','elapsed_s','replacement_prefix_km_s','replacement_core_km_s','replacement_suffix_km_s', ...
 'replacement_visits','replacement_max_miss_m','replacement_min_altitude_km','replacement_days','replacement_burns'});
writetable(tasks,fullfile(folder,'tasks.csv'));control=find(arms=='R');pairs=cell(0,6);
for w=1:W
 for a=find(arms~='R')
  valid=isfinite(cost(w,control))&&isfinite(cost(w,a));gain=NaN;
  if valid,gain=1000*(cost(w,control)-cost(w,a));end
  pairs(end+1,:)={w,d.manifest.source(w),arms(a),valid,available(w,a),gain}; %#ok<AGROW>
 end
end
paired=cell2table(pairs,'VariableNames',{'window','source','arm','complete_pair','replacement_within_budget','net_gain_vs_R_m_s'});
writetable(paired,fullfile(folder,'paired.csv'));report=struct('tasks',tasks,'paired',paired);
if all(ismember('PL',arms))
 ip=find(arms=='P');il=find(arms=='L');lp=cell(W,5);
 for w=1:W
  valid=isfinite(cost(w,ip))&&isfinite(cost(w,il));gain=NaN;
  if valid,gain=1000*(cost(w,ip)-cost(w,il));end
  lp(w,:)={w,valid,available(w,ip),available(w,il),gain};
 end
 report.low_speed_pair=cell2table(lp,'VariableNames',{'window','complete_pair','P_replacement_within_budget','L_replacement_within_budget','gain_vs_P_m_s'});
 writetable(report.low_speed_pair,fullfile(folder,'low_speed_paired.csv'));
 disp(report.low_speed_pair);
end
save(fullfile(folder,'report.mat'),'report');disp(paired);
end
