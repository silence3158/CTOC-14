function report=report_targeted_campaign(label)
%REPORT_TARGETED_CAMPAIGN Report actual outcomes, including failed fragment trials.
root=fileparts(fileparts(mfilename('fullpath')));cd(root);addpath('src');
folder=fullfile(root,'runs','fragments',label);a=load(fullfile(folder,'study.mat'));study=a.study;
a=load(fullfile(folder,'elite.mat'));elite=a.elite;p=a.p;cfg=a.cfg;
v=load(fullfile(folder,'tight_verification.mat'));tight=v.report;assert(tight.passed);
rows=zeros(numel(study.results),12);accepted=zeros(0,6);diagFailures=0;totalEvaluations=0;
for j=1:numel(study.results)
 a=load(fullfile(folder,sprintf('task%02d',j),'result.mat'));o=a.out;h=o.history;
 merge=h(:,3)==2;repair=h(:,3)==1;
 previous=[o.best.independent.total_dv_km_s;h(1:end-1,5)];
 wave=load(fullfile(folder,sprintf('wave%d_source.mat',ceil(j/4))),'source');previous(1)=wave.source.independent.total_dv_km_s;
 failedChecks=isfinite(h(:,4))&h(:,4)<previous-1e-7&~logical(h(:,7));
 rows(j,:)=[j o.config.seed o.elapsed_s size(h,1) sum(merge) sum(merge&isfinite(h(:,4))) ...
  sum(merge&logical(h(:,6))) sum(repair&logical(h(:,6))) sum(failedChecks) h(end,5) h(end,8) ceil(j/4)];
 keep=find(h(:,6)>0);
 accepted=[accepted;[repmat(j,numel(keep),1) h(keep,[2 3 5 8]) (previous(keep)-h(keep,5))*1000]]; %#ok<AGROW>
 for k=1:numel(o.diagnostics)
  d=o.diagnostics{k};
  if isfield(d,'failure_reason'),diagFailures=diagFailures+1;
  else,totalEvaluations=totalEvaluations+d.solver.evaluations;end
 end
end
taskTable=array2table(rows,'VariableNames',{'task','seed','seconds','trials','merge_trials','merge_feasible', ...
 'merge_accepted','connection_accepted','independent_rejected','best_dv_km_s','burns','wave'});
writetable(taskTable,fullfile(folder,'task_summary.csv'));
writetable(array2table(accepted,'VariableNames',{'task','iteration','operator','dv_km_s','burns','saving_m_s'}),fullfile(folder,'accepted.csv'));
s=elite.schedule;ir=elite.independent;
pl=ctocscreen.arcPlan(s,elite.evaluation);ends=cumsum(pl.counts);multi=find(pl.counts>1);
fragmentTable=table('Size',[numel(multi) 5],'VariableTypes',{'double','string','double','double','double'}, ...
 'VariableNames',{'arc','targets','first_visit_s','last_visit_s','max_tight_distance_m'});
for k=1:numel(multi)
 m=multi(k);ix=(ends(m)-pl.counts(m)+1):ends(m);
 fragmentTable(k,:)={m,strjoin(string(pl.ids(ix)),' -> '),pl.times(ix(1)),pl.times(ix(end)),max(tight.event_distances_km(ix))*1000};
end
writetable(fragmentTable,fullfile(folder,'multi_target_arcs.csv'));
writetable(table(s.event_target_ids,s.event_times_s,ir.event_distances_km,tight.event_distances_km, ...
 'VariableNames',{'target','time_s','independent_distance_km','tight_distance_km'}),fullfile(folder,'events.csv'));
writetable(array2table([s.maneuver_times_s s.delta_v_km_s vecnorm(s.delta_v_km_s,2,2)], ...
 'VariableNames',{'time_s','dvx','dvy','dvz','dv_km_s'}),fullfile(folder,'maneuvers.csv'));
[~,winningTask]=min(rows(:,10));elite.winning_task=winningTask;elite.search_config=study.results{winningTask}.config;
save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
report=struct('baseline_dv_km_s',study.baseline,'best_dv_km_s',ir.total_dv_km_s, ...
 'saving_m_s',(study.baseline-ir.total_dv_km_s)*1000,'winning_task',winningTask, ...
 'task_table',taskTable,'total_search_seconds',sum(rows(:,3)),'total_evaluations',totalEvaluations, ...
 'uncaught_trial_failures',diagFailures,'multi_target_arcs',fragmentTable,'tight',tight,'validation_level','two_body_independent_verified');
save(fullfile(folder,'report.mat'),'report');
f=figure('Visible','off','Color','w');tiledlayout(1,2);
nexttile;bar(rows(:,1),(study.baseline-rows(:,10))*1000);xlabel('Task');ylabel('Saving from baseline (m/s)');grid on;
nexttile;plot(0,study.baseline,'ko');hold on;best=study.baseline;yy=best;
for w=1:max(rows(:,12)),best=min(best,min(rows(rows(:,12)==w,10)));yy(end+1)=best;end %#ok<AGROW>
plot(0:max(rows(:,12)),yy,'-o','LineWidth',1.5);xlabel('Completed wave');ylabel('Best verified Delta V (km/s)');grid on;
exportgraphics(f,fullfile(folder,'campaign.png'),'Resolution',160);close(f);
fid=fopen(fullfile(root,'docs','V2_TARGETED_RESULTS_20260921.md'),'w','n','UTF-8');clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# 定向片段搜索批次实测\n\n归档：`runs/fragments/%s`。二体离线筛选，未进行 J2 或 ATK 验证。\n\n',label);
fprintf(fid,'## 完整任务结果\n\n原基线 **%.12f km/s**；本轮最好 **%.12f km/s**，降低 **%.6f m/s（%.4f%%）**。\n\n',study.baseline,ir.total_dv_km_s,report.saving_m_s,100*(1-ir.total_dv_km_s/study.baseline));
fprintf(fid,'35 个不同目标均有实际距离不超过 1 km 的独立积分见证。任务时长 %.9f 天，%d 次脉冲；更严积分全过程最低高度 %.9f km，最大访问距离 %.9f m。\n\n',s.duration_s/86400,numel(s.maneuver_times_s),tight.min_altitude_km,tight.max_distance_km*1000);
fprintf(fid,'初次独立积分最大访问距离 %.9f m；更严积分与初次复核的最大访问距离变化 %.9f m。后者使用 ode113、RelTol=2.3e-14、AbsTol=1e-15、MaxStep=600 s，目标分别从历元积分，固定原脉冲连续重放，不重新瞄准。高度使用有限二体圆锥弧解析极值。\n\n',max(ir.event_distances_km)*1000,tight.max_distance_change_km*1000);
fprintf(fid,'最终保留的多目标自然弧：\n\n');
for k=1:height(fragmentTable)
 fprintf(fid,'- 弧 %d：Target %s；两端访问时刻 %.6f s 至 %.6f s；这些访问之间没有脉冲，更严复核最大距离 %.9f m。\n', ...
  fragmentTable.arc(k),char(fragmentTable.targets(k)),fragmentTable.first_visit_s(k),fragmentTable.last_visit_s(k),fragmentTable.max_tight_distance_m(k));
end
fprintf(fid,'\n');
fprintf(fid,'## 算子与投入\n\n%d 个任务，每批最多 4 个进程、每任务 480 秒软搜索预算。实际合计搜索 %.1f 秒（%.2f 进程分钟），共 %d 次试验、%d 次内层轨迹评价。启动、进程池与复核另计；不把合计进程时间称为墙钟时间。\n\n',height(taskTable),sum(rows(:,3)),sum(rows(:,3))/60,sum(rows(:,4)),totalEvaluations);
fprintf(fid,'合并/重排试验 %d 次，名义全任务通过 %d 次，严格改善并通过独立复核 %d 次；连接窗口优化接受 %d 次。改善候选独立复核拒绝 %d 次，试验级异常 %d 次。\n\n',sum(rows(:,5)),sum(rows(:,6)),sum(rows(:,7)),sum(rows(:,8)),sum(rows(:,9)),diagFailures);
if sum(rows(:,7))==0
 fprintf(fid,'**本轮没有得到可保留的降本多目标合并片段。收益来自连接窗口的时间联合优化，不能将其归因于“一石多鸟”。**\n\n');
else
 fprintf(fid,'本轮存在通过独立复核的合并片段改善。它们的目标组合、访问顺序、实际脉冲和继承关系须结合各任务 result.mat、accepted.csv 及 wave*_source.mat 追溯；不能只以最终脉冲数量推断全部算子贡献。\n\n');
end
fprintf(fid,'|任务|总 ΔV km/s|脉冲数|试验数|合并接受|连接接受|实际秒数|\n|---:|---:|---:|---:|---:|---:|---:|\n');
for j=1:size(rows,1),fprintf(fid,'|%d|%.9f|%d|%d|%d|%d|%.1f|\n',rows(j,[1 10 11 4 7 8 3]));end
fprintf(fid,'\n## 范围及复现\n\n算法与限制见 [V2_TARGETED_SEARCH.md](V2_TARGETED_SEARCH.md)。本轮固定既有自由初轨解，优化目标组合及连接窗口，未穷尽全部拓扑或所有时间范围。旧精英保持原路径不变。任务种子、半径、源码哈希、MAT 检查点、随机状态、候选池与求解器诊断保存在每个 task 目录；源码快照在 source_snapshot。\n\n18 项相关回归测试通过，记录为 `runs/fragments/targeted_tests_20260921.mat`。45 秒烟雾试验记录为 `runs/fragments/targeted_smoke2_20260921`。\n\n最终入口为 `elite.mat`（elite、p、cfg）；CSV 为查看用派生数据，以 MAT 双精度为执行依据。额外验收为 tight_verification.mat，统计为 task_summary.csv、accepted.csv、report.mat，图为 campaign.png。本轮结果不代表全局最优，不直接等同排行榜正式成绩。\n');
disp(report);
end
