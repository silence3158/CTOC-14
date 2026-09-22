function report=report_joint_control_experiment(folder)
%REPORT_JOINT_CONTROL_EXPERIMENT Include missing arms, pair only complete arms.
d=load(fullfile(folder,'design.mat'));groups=d.cfg.groups;K=numel(d.seeds);
rows=cell(K*numel(groups),14);cost=nan(K,numel(groups));at=0;
for k=1:K
 for j=1:numel(groups)
  at=at+1;g=groups(j);file=fullfile(folder,sprintf('seed%02d_%s',k,g),'result.mat');
  status='not_completed';dv=NaN;budgetDv=NaN;elapsed=NaN;evals=NaN;rejects=NaN;qchange=NaN;visits=NaN;alt=NaN;miss=NaN;
  initial=d.seeds{k}.independent.total_dv_km_s;
  if isfile(file)
   a=load(file,'out');o=a.out;status=o.status;elapsed=o.elapsed_s;evals=o.evaluations;rejects=o.validation_rejects;
   ir=o.best.independent;
   if ir.passed
    dv=ir.total_dv_km_s;visits=ir.unique_visit_count;alt=ir.min_altitude_km;miss=max(ir.event_distances_km)*1000;
    dq=o.best.schedule.initial_q-d.seeds{k}.schedule.initial_q;
    dq(5:6)=mod(dq(5:6)+pi,2*pi)-pi;qchange=norm(dq./[10 .001 .001 pi 2*pi 2*pi]);
    budgetDv=o.budget_best.independent.total_dv_km_s;
    if strcmp(status,'completed'),cost(k,j)=budgetDv;end
   end
  end
  rows(at,:)={k,g,status,initial,dv,budgetDv,1000*(initial-budgetDv),elapsed,evals,rejects,qchange,visits,alt,miss};
 end
end
tasks=cell2table(rows,'VariableNames',{'seed','group','status','initial_dv_km_s','best_dv_km_s','budget_dv_km_s', ...
 'gain_m_s','elapsed_s','continuous_and_window_evaluations','continuous_validation_rejects', ...
 'scaled_q_change','verified_visits','min_altitude_km','max_miss_m'});
writetable(tasks,fullfile(folder,'tasks.csv'));
summary=cell(numel(groups),6);
for j=1:numel(groups)
 valid=isfinite(cost(:,j));initial=cellfun(@(s)s.independent.total_dv_km_s,d.seeds)';
 values=cost(valid,j);gains=1000*(initial(valid)-values);
 if isempty(values),best=NaN;med=NaN;gain=NaN;else,best=min(values);med=median(values);gain=median(gains);end
 summary(j,:)={groups(j),sum(valid),K,best,med,gain};
end
summary=cell2table(summary,'VariableNames',{'group','completed','planned','best_km_s','median_km_s','median_gain_m_s'});
writetable(summary,fullfile(folder,'groups.csv'));
pairs=cell(0,7);
for j=2:numel(groups)
 valid=isfinite(cost(:,j-1))&isfinite(cost(:,j));gain=1000*(cost(valid,j-1)-cost(valid,j));
 if isempty(gain),med=NaN;lo=NaN;hi=NaN;else,med=median(gain);lo=min(gain);hi=max(gain);end
 pairs(end+1,:)={sprintf('%c_minus_%c',groups(j-1),groups(j)),sum(valid),med,sum(gain>1),sum(gain<-1),lo,hi}; %#ok<AGROW>
end
pairs=cell2table(pairs,'VariableNames',{'comparison','complete_pairs','median_gain_m_s','wins_over_1m_s','losses_over_1m_s','min_gain_m_s','max_gain_m_s'});
writetable(pairs,fullfile(folder,'paired.csv'));
report=struct('tasks',tasks,'groups',summary,'paired',pairs);
save(fullfile(folder,'report.mat'),'report');disp(summary);disp(pairs);
fig=figure('Visible','off','Position',[100 100 1000 260*ceil(K/2)]);
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
tiledlayout(ceil(K/2),min(2,K));
for k=1:K
 nexttile;hold on;
 for j=1:numel(groups)
  file=fullfile(folder,sprintf('seed%02d_%s',k,groups(j)),'result.mat');
  if isfile(file)
   a=load(file,'out');h=a.out.history;
   stairs(h(:,1),h(:,2),'DisplayName',groups(j),'LineWidth',1.2);
  end
 end
 xline(d.cfg.max_wall_s,'--','Budget','HandleVisibility','off');
 title(sprintf('Paired seed %d',k));xlabel('Elapsed seconds');ylabel('Verified delta-V (km/s)');grid on;
 if ~isempty(findobj(gca,'Type','Stair')),legend('Location','best');end
end
exportgraphics(fig,fullfile(folder,'convergence.png'),'Resolution',140);
end
