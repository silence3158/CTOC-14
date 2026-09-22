function report=report_diverse_campaign(label)
%REPORT_DIVERSE_CAMPAIGN Separate basin exploration from improvement of the incumbent.
root=fileparts(fileparts(mfilename('fullpath')));cd(root);addpath('src');
folder=fullfile(root,'runs','fragments',label);a=load(fullfile(folder,'study.mat'));study=a.study;
a=load(fullfile(folder,'elite.mat'));elite=a.elite;p=a.p;cfg=a.cfg;
seeds=load(fullfile(folder,'seeds.mat'));v=load(fullfile(folder,'tight_verification.mat'));tight=v.report;assert(tight.passed);
rows=zeros(numel(study.results),11);allOrders=seeds.orders;exceptions=0;totalEvaluations=0;
for j=1:numel(study.results)
 a=load(fullfile(folder,sprintf('task%02d',j),'result.mat'));o=a.out;h=o.history;
 restarts=sum(cellfun(@(r)isfield(r,'order_hamming'),o.restarts));
 rows(j,:)=[j seeds.sources{j}.independent.total_dv_km_s o.best.independent.total_dv_km_s ...
 o.elapsed_s size(h,1) sum(h(:,3)==2) sum(h(:,3)==2&h(:,6)>0) sum(h(:,3)==1&h(:,6)>0) ...
 restarts numel(o.best.schedule.maneuver_times_s) o.config.seed];
 allOrders(end+1,:)=o.current.schedule.event_target_ids'; %#ok<AGROW>
 for k=1:numel(o.diagnostics)
  d=o.diagnostics{k};if isfield(d,'failure_reason'),exceptions=exceptions+1;
  else,totalEvaluations=totalEvaluations+d.solver.evaluations;end
 end
end
taskTable=array2table(rows,'VariableNames',{'task','initial_dv_km_s','best_dv_km_s','seconds','trials', ...
 'merge_trials','merge_accepted_in_chain','connection_accepted_in_chain','successful_restarts','best_burns','seed'});
writetable(taskTable,fullfile(folder,'task_summary.csv'));
s=elite.schedule;ir=elite.independent;
writetable(table(s.event_target_ids,s.event_times_s,ir.event_distances_km,tight.event_distances_km, ...
 'VariableNames',{'target','time_s','independent_distance_km','tight_distance_km'}),fullfile(folder,'events.csv'));
writetable(array2table([s.maneuver_times_s s.delta_v_km_s vecnorm(s.delta_v_km_s,2,2)], ...
 'VariableNames',{'time_s','dvx','dvy','dvz','dv_km_s'}),fullfile(folder,'maneuvers.csv'));
pl=ctocscreen.arcPlan(s,elite.evaluation);ends=cumsum(pl.counts);multi=find(pl.counts>1);
fragmentTable=table('Size',[numel(multi) 5],'VariableTypes',{'double','string','double','double','double'}, ...
 'VariableNames',{'arc','targets','first_visit_s','last_visit_s','max_tight_distance_m'});
for k=1:numel(multi)
 m=multi(k);ix=(ends(m)-pl.counts(m)+1):ends(m);
 fragmentTable(k,:)={m,strjoin(string(pl.ids(ix)),' -> '),pl.times(ix(1)),pl.times(ix(end)),max(tight.event_distances_km(ix))*1000};
end
writetable(fragmentTable,fullfile(folder,'multi_target_arcs.csv'));
report=struct('baseline_dv_km_s',study.baseline,'best_dv_km_s',ir.total_dv_km_s, ...
 'saving_m_s',(study.baseline-ir.total_dv_km_s)*1000,'task_table',taskTable, ...
 'distinct_initial_orders',size(unique(seeds.orders,'rows'),1),'initial_and_final_orders',size(unique(allOrders,'rows'),1), ...
 'successful_restarts',sum(rows(:,9)),'total_search_seconds',sum(rows(:,4)), ...
 'total_evaluations',totalEvaluations,'trial_exceptions',exceptions,'multi_target_arcs',fragmentTable,'tight',tight);
save(fullfile(folder,'report.mat'),'report');
f=figure('Visible','off','Color','w');tiledlayout(1,2);
nexttile;plot(rows(:,1),rows(:,2),'o--',rows(:,1),rows(:,3),'o-');hold on;yline(study.baseline,'k:');
xlabel('Task');ylabel('Delta V (km/s)');legend('Initial seed','Task best','Retained incumbent','Location','best');grid on;
nexttile;bar(rows(:,1),rows(:,9));xlabel('Task');ylabel('Successful basin restarts');grid on;
exportgraphics(f,fullfile(folder,'campaign.png'),'Resolution',160);close(f);
fid=fopen(fullfile(root,'docs','V2_DIVERSE_RESULTS_20260921.md'),'w','n','UTF-8');clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# 多起点与跳出局部区域：实测记录\n\n实际批次：`runs/fragments/%s`。全部结论限于二体，未运行 J2 或 ATK。\n\n',label);
fprintf(fid,'## 用户纠正后的算法改变\n\n原批次只换 RNG 种子、各波共享最好解，主要在同一局部区域下降。按用户实时指示，`targeted16_20260921` 在 4 个任务完成、4 个运行途中停止，剩余 8 个未启动。已保存检查点，不称其为完成的 16 任务。该阶段从 14.063572065649 降至 13.976837787874 km/s，得到 Target20 → Target33 双目标自然弧；大部分降幅来自连接时刻优化，接受合并那一步的净收益约 1.8–2.9 m/s。\n\n');
fprintf(fid,'新批次先独立验证 %d 个访问顺序两两不同的完整任务，再按 4 进程、4 波、每任务 480 秒软预算优化。种子来源和随机状态在 seeds.mat；实际顺序差异及时间 RMS 在 seed_manifest.csv。各波不重置到共同精英。工作解可暂时更贵，历史最好解独立保留；停滞或周期触发换序、块搬移/反转、时间重分配、重接分支与拆开旧片段。\n\n',report.distinct_initial_orders);
fprintf(fid,'原解的长滑行和分支参考与原访问顺序绑定，不直接移植到强扰动后。构造新种子时先使用访问后立即连接、按实际到达速度选择局部低脉冲分支；局部优化时重新开放滑行时间。它是初始化方法，不是正式模型约束，也不意味着逐段贪心为全局最优。\n\n');
fprintf(fid,'## 实际结果\n\n新批次进入时保底 **%.12f km/s**，最终保留 **%.12f km/s**，本轮新降低 **%.6f m/s**。\n\n',study.baseline,ir.total_dv_km_s,report.saving_m_s);
if report.saving_m_s<=1e-4
 fprintf(fid,'**本轮确实探索了不同访问顺序并执行了上坡重启，但尚未超过保底解。不能把“已经跳出原区域”写成“已经找到更好的区域”。**\n\n');
end
fprintf(fid,'实际成功重启 %d 次，合计搜索 %.1f 秒（%.2f 进程分钟），试验 %d 次、内层轨迹评价 %d 次。初始与各任务最终工作解至少覆盖 %d 个不同顺序（此数不包含所有中途访问过的顺序）。试验级异常 %d 次。种子构造、MATLAB 启动、进程池与最终复核不计入单任务搜索时间。\n\n',report.successful_restarts,report.total_search_seconds,report.total_search_seconds/60,sum(rows(:,5)),totalEvaluations,report.initial_and_final_orders,exceptions);
fprintf(fid,'最终 35 个不同目标有实际独立积分访问见证；%d 次脉冲，时长 %.9f 天。更严复核最大访问距离 %.9f m，最低高度 %.9f km，和初次独立复核的最大距离变化 %.9f m。高度覆盖全部有限二体弧的解析径向极值。更严复核固定原脉冲连续重放，目标分别自历元积分，ode113 的 RelTol=2.3e-14、AbsTol=1e-15、MaxStep=600 s。不声称穷尽首次进入时刻或额外近遇。\n\n',numel(s.maneuver_times_s),s.duration_s/86400,tight.max_distance_km*1000,tight.min_altitude_km,tight.max_distance_change_km*1000);
fprintf(fid,'最终多目标自然弧：\n\n');
for k=1:height(fragmentTable),fprintf(fid,'- Target %s，访问时刻 %.6f s 至 %.6f s，期间没有脉冲。\n',char(fragmentTable.targets(k)),fragmentTable.first_visit_s(k),fragmentTable.last_visit_s(k));end
fprintf(fid,'\n|任务|起始 ΔV|任务最好 ΔV|成功重启|最好脉冲数|实际秒数|\n|---:|---:|---:|---:|---:|---:|\n');
for j=1:size(rows,1),fprintf(fid,'|%d|%.9f|%.9f|%d|%d|%.1f|\n',rows(j,[1 2 3 9 10 4]));end
fprintf(fid,'\n## 证据及边界\n\n`elite.mat` 为可接入结果（elite、p、cfg）；`tight_verification.mat` 为严格复核；`task_summary.csv`、`report.mat`、`campaign.png` 为汇总；每个 task 的 result.mat 保留当前/历史最好解、重启诊断、候选组合、求解器信息、随机状态与配置。CSV 仅便于查看，执行以 MAT 双精度值为准。源码快照在 source_snapshot。18 项回归测试与非局部种子复核记录为 `runs/fragments/diverse_final_tests_20260921.mat`，上坡工作状态与最好解隔离的实测为 `runs/fragments/diverse_escape_smoke_20260921`。\n\n新批次使用不同种子携带的初轨，但每条搜索链内未连续优化初轨六自由度；单次修正的时间范围、搜索圈数、排序规则、重启阈值和种子质量筛选均为算法近似。这里只能报告有限预算的搜索结果，不能证明局部最优或全局最优。\n');
disp(report);
end
