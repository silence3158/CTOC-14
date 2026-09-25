function summary=summarize_zhang_tail27()
%SUMMARIZE_ZHANG_TAIL27 Compare archived independent results; no optimization.
sim=fileparts(fileparts(mfilename('fullpath'))); base=fullfile(sim,'runs/v3/diagnostics');
out=fullfile(base,'zhang_tail27_comparison_03');
assert(~isfolder(out),'Do not overwrite comparison evidence.'); mkdir(out);
a=load(fullfile(base,'tail27_20260925_method_search/comparison.mat'),'result');
b=load(fullfile(base,'zhang_dp_tail27_visit_bound_02/experiment.mat'),'result');
d=load(fullfile(base,'zhang_delayed_tail27_delayed_03/experiment.mat'),'result');
assert(b.result.source_unchanged&&d.result.source_unchanged);
schedules={a.result.reconstructed_schedule,b.result.best_new.schedule,d.result.best_new.schedule};
checks={a.result.reconstructed_verification,b.result.best_new.verification,d.result.best_new.verification};
names={'teammate_reconstructed';'visit_bound_dp';'delayed_finite_beam'};
data=zeros(3,8); order=a.result.plan(:,1);
for k=1:3
 s=schedules{k}; v=checks{k}; raw=sum(vecnorm(s.delta_v_km_s,2,2));
 assert(v.passed&&v.visit_count==35&&v.targets_independently_propagated);
 assert(abs(raw-v.total_dv_km_s)<1e-10);
 assert(all(diff(s.witness_times_s(order))>0));
 waits=s.maneuver_times_s-[0;s.witness_times_s(order(1:end-1))];
 assert(all(waits>=-1e-9));
 data(k,:)=[raw,v.visit_count,max(v.distance_km),v.min_altitude_lower_km, ...
  s.duration_s/86400,numel(s.maneuver_times_s),sum(waits(2:end)>60),sum(waits)/86400];
 burns=array2table([(1:35).',order,s.maneuver_times_s,s.witness_times_s(order),waits, ...
  s.delta_v_km_s,vecnorm(s.delta_v_km_s,2,2)],'VariableNames', ...
  {'leg','target','burn_s','arrival_s','wait_s','dv_x_km_s','dv_y_km_s','dv_z_km_s','dv_km_s'});
 writetable(burns,fullfile(out,[names{k} '_burns.csv']));
end
summary=array2table(data,'VariableNames',{'raw_dv_km_s','visits','max_distance_km', ...
 'min_altitude_lower_km','duration_days','burn_count','post_visit_waits_over_60s','total_wait_days'});
summary=addvars(summary,names,'Before',1,'NewVariableNames','variant');
summary.experiment_elapsed_s=[NaN;b.result.elapsed_s;d.result.elapsed_s];
summary.improvement_km_s=data(1,1)-data(:,1);
writetable(summary,fullfile(out,'summary.csv')); save(fullfile(out,'summary.mat'),'summary');
fig=figure('Visible','off','Color','w','Position',[100 100 1100 460]);
tiledlayout(1,2,'Padding','compact');
nexttile; bar(data(:,1),'FaceColor',[.18 .53 .43]);
xticks(1:3); xticklabels({'Teammate','Visit-bound DP','Delayed beam'}); ylabel('Raw total delta-V (km/s)');
for k=1:3, text(k,data(k,1)+.15,sprintf('%.6f',data(k,1)),'HorizontalAlignment','center'); end
ylim([0,max(data(:,1))+2]); grid on;
nexttile; hold on;
for k=1:3, plot(1:35,cumsum(vecnorm(schedules{k}.delta_v_km_s,2,2)),'LineWidth',1.5); end
xlabel('Planned visit'); ylabel('Cumulative raw delta-V (km/s)');
legend({'Teammate','Visit-bound DP','Delayed beam'},'Location','northwest'); grid on;
exportgraphics(fig,fullfile(out,'comparison.png')); close(fig);
disp(summary); fprintf('COMPARISON %s\n',out);
end
