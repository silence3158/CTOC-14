function report=report_diverse_campaign(label)
%REPORT_DIVERSE_CAMPAIGN Separate basin exploration from improvement of the incumbent.
root=fileparts(fileparts(mfilename('fullpath')));cd(root);addpath('src');
folder=fullfile(root,'runs','fragments',label);a=load(fullfile(folder,'study.mat'));study=a.study;
a=load(fullfile(folder,'elite.mat'));elite=a.elite;p=a.p;cfg=a.cfg;
seeds=load(fullfile(folder,'seeds.mat'));v=load(fullfile(folder,'tight_verification.mat'));tight=v.report;assert(tight.passed);
origin=load(fullfile(folder,'config.mat'),'base');baseTimes=zeros(35,1);baseTimes(origin.base.schedule.event_target_ids)=origin.base.schedule.event_times_s;
seedDetails=zeros(numel(seeds.sources),9);
for j=1:numel(seeds.sources)
 ss=seeds.sources{j}.schedule;targetTimes=zeros(35,1);targetTimes(ss.event_target_ids)=ss.event_times_s;
 seedDetails(j,:)=[j ss.initial_q(:)' sqrt(mean((ss.event_times_s-origin.base.schedule.event_times_s).^2)) sqrt(mean((targetTimes-baseTimes).^2))];
end
writetable(array2table(seedDetails,'VariableNames',{'task','a_km','ec','es','inclination_rad','raan_rad','u_rad','slot_time_rms_s','per_target_time_rms_s'}),fullfile(folder,'seed_details.csv'));
rows=zeros(numel(study.results),11);allOrders=seeds.orders;exceptions=0;totalEvaluations=0;fragments={};
for j=1:numel(study.results)
 a=load(fullfile(folder,sprintf('task%02d',j),'result.mat'));o=a.out;h=o.history;
 restarts=sum(cellfun(@(r)isfield(r,'order_hamming'),o.restarts));
 rows(j,:)=[j seeds.sources{j}.independent.total_dv_km_s o.best.independent.total_dv_km_s ...
 o.elapsed_s numel(o.diagnostics) sum(h(:,3)==2) sum(h(:,3)==2&h(:,6)>0) sum(h(:,3)==1&h(:,6)>0) ...
 restarts numel(o.best.schedule.maneuver_times_s) o.config.seed];
 allOrders(end+1,:)=o.current.schedule.event_target_ids'; %#ok<AGROW>
 fragments=[fragments extractFragments(o.best,j,'task_best',p) extractFragments(o.current,j,'final_working',p)]; %#ok<AGROW>
 for k=1:numel(o.diagnostics)
  d=o.diagnostics{k};if isfield(d,'failure_reason'),exceptions=exceptions+1;
  else,totalEvaluations=totalEvaluations+d.solver.evaluations;end
 end
end
fragments=[extractFragments(elite,0,'retained_global_best',p) fragments];
rawCatalogCount=numel(fragments);distinct={};
for k=1:numel(fragments)
 ff=fragments{k};origin=struct('task',ff.task,'source',ff.source);match=[];
 for j=1:numel(distinct)
  gg=distinct{j};
  if isequal(ff.event_target_ids,gg.event_target_ids)&&isequal(ff.event_times_s,gg.event_times_s)&& ...
    isequal(ff.t_in_s,gg.t_in_s)&&isequal(ff.state_in,gg.state_in)&&isequal(ff.delta_v_km_s,gg.delta_v_km_s),match=j;break;end
 end
 if isempty(match),ff.origins={origin};distinct{end+1}=ff; %#ok<AGROW>
 else,distinct{match}.origins{end+1}=origin;distinct{match}.source=[distinct{match}.source ' + ' ff.source];end
end
fragments=distinct;
save(fullfile(folder,'fragment_catalog.mat'),'fragments','p','-v7.3');
catalog=table('Size',[numel(fragments) 7],'VariableTypes',{'double','string','string','double','double','double','double'}, ...
 'VariableNames',{'task','source','targets','entry_s','exit_s','local_dv_km_s','whole_mission_dv_km_s'});
for k=1:numel(fragments)
 ff=fragments{k};catalog(k,:)={ff.task,string(ff.source),strjoin(string(ff.event_target_ids),' -> '),ff.t_in_s,ff.t_out_s,ff.total_dv_km_s,ff.whole_mission_dv_km_s};
end
writetable(catalog,fullfile(folder,'fragment_catalog.csv'));
taskTable=array2table(rows,'VariableNames',{'task','initial_dv_km_s','best_dv_km_s','seconds','trials', ...
 'merge_trials','merge_accepted_in_chain','connection_accepted_in_chain','successful_restarts','best_burns','seed'});
writetable(taskTable,fullfile(folder,'task_summary.csv'));
s=elite.schedule;ir=elite.independent;
if ir.total_dv_km_s>=study.baseline-1e-7
 elite.provenance=struct('origin','retained prior verified incumbent','source','runs/fragments/targeted16_20260921/elite.mat');
else
 elite.provenance=struct('origin','improved in diverse campaign','task',elite.winning_task);
end
elite.provenance.final_verification_source_hash=ctocscreen.implementationHash();
save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
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
 'total_evaluations',totalEvaluations,'trial_exceptions',exceptions,'multi_target_arcs',fragmentTable, ...
 'raw_catalog_records',rawCatalogCount,'fragment_catalog_records',numel(fragments),'distinct_fragment_combinations',numel(unique(catalog.targets)),'tight',tight);
save(fullfile(folder,'report.mat'),'report');
f=figure('Visible','off','Color','w','Position',[100 100 1150 500]);tiledlayout(1,2);
nexttile;plot(rows(:,1),rows(:,2),'o--',rows(:,1),rows(:,3),'o-');hold on;yline(study.baseline,'k:');
xlabel('Task');ylabel('Delta V (km/s)');legend('Initial seed','Task best','Retained incumbent','Location','northoutside','NumColumns',3);grid on;
nexttile;bar(rows(:,1),rows(:,9));xlabel('Task');ylabel('Successful basin restarts');yticks(0:max(rows(:,9)));grid on;
exportgraphics(f,fullfile(folder,'campaign.png'),'Resolution',160);close(f);
fid=fopen(fullfile(root,'docs','V2_DIVERSE_RESULTS_20260921.md'),'w','n','UTF-8');clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# 多起点与跳出局部区域：实测记录\n\n实际批次：`runs/fragments/%s`。全部结论限于二体，未运行 J2 或 ATK。\n\n',label);
fprintf(fid,'## 用户纠正后的算法改变\n\n原批次只换 RNG 种子、各波共享最好解，主要在同一局部区域下降。按用户实时指示，`targeted16_20260921` 在 4 个任务完成、4 个运行途中停止，剩余 8 个未启动。已保存检查点，不称其为完成的 16 任务。该阶段从 14.063572065649 降至 13.976837787874 km/s，得到 Target20 → Target33 双目标自然弧；大部分降幅来自连接时刻优化，接受合并那一步的净收益约 1.8–2.9 m/s。\n\n');
fprintf(fid,'新批次先独立验证 %d 个访问顺序两两不同的完整任务，再按 4 进程、4 波、每任务 480 秒软预算优化。种子来源和随机状态在 seeds.mat；实际顺序差异及时间 RMS 在 seed_manifest.csv。各波不重置到共同精英。工作解可暂时更贵，历史最好解独立保留；停滞或周期触发换序、块搬移/反转、时间重分配、重接分支与拆开旧片段。\n\n',report.distinct_initial_orders);
fprintf(fid,'原解的长滑行和分支参考与原访问顺序绑定，不直接移植到强扰动后。构造新种子时先使用访问后立即连接、按实际到达速度选择局部低脉冲分支；局部优化时重新开放滑行时间。它是初始化方法，不是正式模型约束，也不意味着逐段贪心为全局最优。\n\n');
fprintf(fid,'历史种子的零等待曾带有约 1e-11 s 的负舍入残差，导致个别早期试验无法重建。第一波完成后、第二波运行途中有序暂停，在既有同刻事件容差 1e-7 s 内将这种负零规范为 0。通过新增回归和全部 16 种子重建检查后，已完成的 4 个任务不重算，其余任务从检查点工作解、最好解与 RNG 状态继续，仅补足每任务总计 480 s 预算。原始输入和原始脉冲归档未改写；规范后的候选仍须独立积分。任务 5–8 的 phase1_result.mat 与 *_continued 目录保留分阶段证据，合并结果记录两版源码哈希；source_snapshot_time_fix 保留修复版本。\n\n');
fprintf(fid,'## 实际结果\n\n新批次进入时保底 **%.12f km/s**，最终保留 **%.12f km/s**，本轮新降低 **%.6f m/s**。\n\n',study.baseline,ir.total_dv_km_s,report.saving_m_s);
if report.saving_m_s<=1e-4
 fprintf(fid,'**本轮确实探索了不同访问顺序并执行了上坡重启，但尚未超过保底解。不能把“已经跳出原区域”写成“已经找到更好的区域”。**\n\n');
end
fprintf(fid,'实际成功重启 %d 次，合计搜索 %.1f 秒（%.2f 进程分钟），共 %d 条试验诊断，记录到内层轨迹评价 %d 次。初始与各任务最终工作解至少覆盖 %d 个不同顺序（此数不包含所有中途访问过的顺序）。其中有 %d 次试验级异常：均为任务 9 的候选未能同时通过独立积分与快速固定脉冲重放，按保守策略拒绝，未覆盖精英；未继续诊断两种传播结果差异的具体原因。种子构造、MATLAB 启动、进程池与最终复核不计入单任务搜索时间。\n\n',report.successful_restarts,report.total_search_seconds,report.total_search_seconds/60,sum(rows(:,5)),totalEvaluations,report.initial_and_final_orders,exceptions);
fprintf(fid,'最终 35 个不同目标有实际独立积分访问见证；%d 次脉冲，时长 %.9f 天。更严复核最大访问距离 %.9f m，最低高度 %.9f km，和初次独立复核的最大距离变化 %.9f m。高度覆盖全部有限二体弧的解析径向极值。更严复核固定原脉冲连续重放，目标分别自历元积分，ode113 的 RelTol=2.3e-14、AbsTol=1e-15、MaxStep=600 s。不声称穷尽首次进入时刻或额外近遇。\n\n',numel(s.maneuver_times_s),s.duration_s/86400,tight.max_distance_km*1000,tight.min_altitude_km,tight.max_distance_change_km*1000);
fprintf(fid,'最终多目标自然弧：\n\n');
for k=1:height(fragmentTable),fprintf(fid,'- Target %s，访问时刻 %.6f s 至 %.6f s，期间没有脉冲。\n',char(fragmentTable.targets(k)),fragmentTable.first_visit_s(k),fragmentTable.last_visit_s(k));end
fprintf(fid,'\n另外，从各任务最好解与最终工作解提取了 %d 条原始多目标片段记录，合并完全相同的轨迹后保留 %d 条记录、%d 种不同目标组合，保存为 fragment_catalog.mat / fragment_catalog.csv。每条 MAT 记录保留真实入口/出口状态、时刻、已访问掩码、原脉冲、实际访问见证及其所在完整任务的成本；origins 记录合并来源。相同目标组合的不同边界状态分别保留；不能脱离这些条件只按目标编号拼接，也不能把局部脉冲成本当成加入另一个任务后的净收益。\n',rawCatalogCount,numel(fragments),numel(unique(catalog.targets)));
fprintf(fid,'\n|任务|起始 ΔV|任务最好 ΔV|成功重启|最好脉冲数|实际秒数|\n|---:|---:|---:|---:|---:|---:|\n');
for j=1:size(rows,1),fprintf(fid,'|%d|%.9f|%.9f|%d|%d|%.1f|\n',rows(j,[1 2 3 9 10 4]));end
fprintf(fid,'\n## 证据及边界\n\n`elite.mat` 为可接入结果（elite、p、cfg）；`tight_verification.mat` 为严格复核；`task_summary.csv`、`report.mat`、`campaign.png` 为汇总；每个 task 的 result.mat 保留当前/历史最好解、重启诊断、候选组合、求解器信息、随机状态与配置。CSV 仅便于查看，执行以 MAT 双精度值为准。源码快照在 source_snapshot 和 source_snapshot_time_fix。初版 18 项回归与非局部种子复核记录为 `runs/fragments/diverse_final_tests_20260921.mat`；修复后的 19 项回归及全部 16 种子重建检查记录为 `runs/fragments/diverse_time_fix_tests_20260921.mat`。上坡工作状态与最好解隔离的实测为 `runs/fragments/diverse_escape_smoke_20260921`。\n\n新批次使用不同种子携带的初轨，但每条搜索链内未连续优化初轨六自由度；单次修正的时间范围、搜索圈数、排序规则、重启阈值和种子质量筛选均为算法近似。这里只能报告有限预算的搜索结果，不能证明局部最优或全局最优。\n');
disp(report);
end
function entries=extractFragments(source,task,kind,p)
entries={};assert(source.independent.passed);s=source.schedule;ir=source.independent;
pl=ctocscreen.arcPlan(s,source.evaluation);ends=cumsum(pl.counts);
for m=find(pl.counts>1)'
 ix=(ends(m)-pl.counts(m)+1):ends(m);mask=false(1,35);mask(s.event_target_ids(s.event_times_s<=s.maneuver_times_s(m)))=true;
 z=[0;diff([s.maneuver_times_s(m);s.event_times_s(ix)])]';z=[z s.delta_v_km_s(m,:)];
 local=ctocscreen.replayMultiFragment(ir.preburn_states(m,:),s.maneuver_times_s(m),s.event_target_ids(ix),pl.counts(m),z,p);
 assert(all(local.event_distances_km<=1)&&local.min_altitude_km>=200);
 entries{end+1}=struct('task',task,'source',kind,'schema_version','fragment_record_v2','dynamics_id','two_body', ...
  'counts',pl.counts(m),'z',z,'local_replay_distances_km',local.event_distances_km, ...
  't_in_s',s.maneuver_times_s(m),'state_in',ir.preburn_states(m,:),'visited_mask',mask, ...
  'maneuver_times_s',s.maneuver_times_s(m),'delta_v_km_s',s.delta_v_km_s(m,:), ...
  'event_target_ids',s.event_target_ids(ix),'event_times_s',s.event_times_s(ix),'event_states',ir.event_states(ix,:), ...
  'event_distances_km',ir.event_distances_km(ix),'t_out_s',s.event_times_s(ix(end)), ...
  'state_out',ir.event_states(ix(end),:),'total_dv_km_s',norm(s.delta_v_km_s(m,:)), ...
  'whole_mission_dv_km_s',ir.total_dv_km_s,'validation_level','part_of_fixed_impulse_independently_verified_mission'); %#ok<AGROW>
end
end
