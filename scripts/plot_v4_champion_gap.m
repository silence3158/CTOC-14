function plot_v4_champion_gap
%PLOT_V4_CHAMPION_GAP Plot archived fixed-control comparison, no integration.
sim=fileparts(fileparts(mfilename('fullpath')));
data=load(fullfile(sim,'runs/v4/diagnostics/champion_gap_20260928_02/gap_analysis.mat'));
r=data.report; folder=fullfile(sim,'docs/assets'); if ~isfolder(folder), mkdir(folder); end
labels={'冠军（缓存状态反推）','第二名（原始脉冲）','队友解（本地重建）','本方冷启动 30 分钟'};
colors=[.12 .48 .31;.12 .4 .72;.68 .45 .1;.75 .22 .22];
fig=figure('Visible','off','Color','w','Position',[80 80 1200 850]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for k=1:4
 d=r.trajectories{k}; plot(d.visits.order,d.visits.cumulative_dv,'LineWidth',2,'Color',colors(k,:));
end
grid on; xlim([1 35]); ylabel('累计 ΔV（km/s）'); xlabel('实际访问次序');
title('差距随任务推进而扩大：冠军最后 5 个目标仅再花约 0.305 km/s');
legend(labels,'Location','northwest','FontSize',10); ylim([0 19]);
nexttile; hold on; d=r.trajectories{1}; b=d.burns; major=b.dv_km_s>.001;
stem(b.day(major),b.dv_km_s(major),'Color',colors(1,:),'Marker','none','LineWidth',1.3);
scatter(d.visits.day,repmat(-.045,35,1),22,[.16 .3 .65],'filled');
plot([0 10],[0 0],'Color',[.65 .65 .65]); xlim([0 10]); ylim([-.13 1.7]); grid on;
ylabel('单次 ΔV（km/s）'); xlabel('任务时间（天）');
title('冠军：40 次 >1 m/s 机动中，31 次距最近访问超过 60 秒');
legend({'较大机动','35 次目标访问'},'Location','northeast');
nexttile; hold on; d=r.trajectories{4}; b=d.burns; major=b.dv_km_s>.001;
stem(b.day(major),b.dv_km_s(major),'Color',colors(4,:),'Marker','none','LineWidth',1.3);
scatter(d.visits.day,repmat(-.07,35,1),22,[.16 .3 .65],'filled');
plot([0 10],[0 0],'Color',[.65 .65 .65]); xlim([0 10]); ylim([-.18 3.3]); grid on;
ylabel('单次 ΔV（km/s）'); xlabel('任务时间（天）');
title('本方：34 次机动中仅 9 次距访问超过 60 秒，最后 5 个目标再花约 5.757 km/s');
legend({'较大机动','35 次目标访问'},'Location','northwest');
set(findall(fig,'-property','FontName'),'FontName','Microsoft YaHei');
set(findall(fig,'Type','axes'),'FontSize',11);
exportgraphics(fig,fullfile(folder,'V4_CHAMPION_GAP_20260928.png'),'Resolution',160);
close(fig);
end
