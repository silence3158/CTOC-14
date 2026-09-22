function report_fragment_pilot(label)
%REPORT_FRAGMENT_PILOT Reproducible tabular summary and scientific figures.
root=fileparts(fileparts(mfilename('fullpath'))); cd(root); addpath('src');
folder=fullfile('runs/fragments',label); a=load(fullfile(folder,'batch.mat')); batch=a.batch;
n=numel(batch.results); rows=zeros(n,14);
fragment_archive={}; screened_archive={}; export_archive={};
for j=1:n
 o=batch.results{j}; st=o.statistics;
 fragment_archive=[fragment_archive,o.fragment_archive]; %#ok<AGROW>
 screened_archive{end+1}=struct('schedule',o.schedule,'evaluation',o.evaluation); %#ok<AGROW>
 if o.independent.passed
  export_archive{end+1}=struct('schedule',o.schedule,'evaluation',o.independent); %#ok<AGROW>
 end
 rows(j,:)=[j,o.config.seed,o.search_elapsed_s,o.evaluation.total_dv_km_s, ...
  o.independent.passed,st.shots,st.gated,st.fragments,st.reconnected, ...
  st.s2_improvements,st.s1_improvements,st.solver_evaluations, ...
  size(o.schedule.delta_v_km_s,1),sum(o.plan.waits(2:end)>1)];
end
names={'job','seed','search_s','dv_km_s','verified','shots','gated','S2_feasible', ...
 'reconnected','S2_improvements','S1_improvements','solver_evaluations','burn_count','delayed_burns'};
tab=array2table(rows,'VariableNames',names); writetable(tab,fullfile(folder,'summary.csv')); disp(tab);
save(fullfile(folder,'fragment_archive.mat'),'fragment_archive','-v7.3');
save(fullfile(folder,'screened_archive.mat'),'screened_archive','-v7.3');
save(fullfile(folder,'export_archive.mat'),'export_archive','-v7.3');
f=figure('Visible','off','Color','w'); hold on;
for j=1:n
 o=batch.results{j}; plot(o.history(:,1),o.history(:,2),'DisplayName',sprintf('seed %d',o.config.seed));
end
yline(batch.baseline_dv_km_s,'--','V1 baseline','DisplayName','V1 baseline'); xlabel('Search wall time (s)'); ylabel('Screened total delta-V (km/s)');
legend('Location','best'); grid on; exportgraphics(f,fullfile(folder,'convergence.png'),'Resolution',160); close(f);
if batch.best_job==0, return; end
o=batch.results{batch.best_job}; s=o.schedule; r=o.independent;
saved=load(fullfile(folder,'elite.mat'),'p'); p=saved.p;
writetable(table(s.maneuver_times_s,s.delta_v_km_s(:,1),s.delta_v_km_s(:,2),s.delta_v_km_s(:,3), ...
 vecnorm(s.delta_v_km_s,2,2),'VariableNames',{'time_s','dvx','dvy','dvz','dv_km_s'}),fullfile(folder,'maneuvers.csv'));
writetable(table(s.event_target_ids,s.event_times_s,r.event_distances_km, ...
 'VariableNames',{'target','witness_time_s','distance_km'}),fullfile(folder,'events.csv'));
f=figure('Visible','off','Color','w');
tiledlayout(2,1); nexttile; stem(s.maneuver_times_s/86400,vecnorm(s.delta_v_km_s,2,2),'filled');
xlabel('Mission time (days)'); ylabel('Impulse (km/s)'); grid on;
nexttile; scatter(s.event_times_s/86400,s.event_target_ids,25,'filled');
xlabel('Mission time (days)'); ylabel('Visited target ID'); grid on;
exportgraphics(f,fullfile(folder,'events_and_burns.png'),'Resolution',160); close(f);
tt=[]; hh=[];
for k=1:size(r.arc_times_s,1)
 times=linspace(r.arc_times_s(k,1),r.arc_times_s(k,2),40);
 heights=zeros(size(times));
 for j=1:numel(times)
  state=ctocscreen.propagateTwoBody(r.arc_states(k,:),times(j)-times(1),p.mu_km3_s2);
  heights(j)=norm(state(1:3))-p.re_km;
 end
 tt=[tt times]; hh=[hh heights]; %#ok<AGROW>
end
f=figure('Visible','off','Color','w'); plot(tt/86400,hh); hold on; yline(200,'--');
xlabel('Mission time (days)'); ylabel('Altitude (km)'); grid on;
title(sprintf('Analytic arc minimum: %.3f km',r.min_altitude_km));
exportgraphics(f,fullfile(folder,'altitude.png'),'Resolution',160); close(f);
pair=find(o.plan.skip,1);
if ~isempty(pair)
 eventTimes=s.event_times_s(pair-1:pair); ids=s.event_target_ids(pair-1:pair);
 burn=find(s.maneuver_times_s<eventTimes(1),1,'last');
 state=r.preburn_states(burn,:); state(4:6)=state(4:6)+s.delta_v_km_s(burn,:);
 start=s.maneuver_times_s(burn); times=sort(unique([linspace(start,eventTimes(2),500) eventTimes']));
 distances=zeros(numel(times),2);
 for j=1:numel(times)
  xx=ctocscreen.propagateTwoBody(state,times(j)-start,p.mu_km3_s2);
  targets=ctocscreen.targetStates(p,ids,times(j)); distances(j,:)=vecnorm(targets(:,1:3)-xx(1:3),2,2)';
 end
 f=figure('Visible','off','Color','w'); semilogy((times-start)/3600,max(distances,1e-6)); hold on; yline(1,'--');
 xlabel('Time since single impulse (hours)'); ylabel('Distance to target (km)');
 legend(compose('Target %d',ids),'Location','best'); grid on;
 exportgraphics(f,fullfile(folder,'double_flyby.png'),'Resolution',160); close(f);
end
fprintf('BEST %.12f gain%.6f%% visits%d max_error_m%.8f min_alt%.6f duration_days%.9f burns%d\n', ...
 r.total_dv_km_s,100*(1-r.total_dv_km_s/batch.baseline_dv_km_s),r.unique_visit_count, ...
 max(r.event_distances_km)*1000,r.min_altitude_km,s.duration_s/86400,size(s.delta_v_km_s,1));
fid=fopen(fullfile(folder,'summary.md'),'w','n','UTF-8'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# V2 pilot results\n\nScope: %s\n\n',batch.config.scope);
fprintf(fid,'Baseline: %.12f km/s. Best independently verified: %.12f km/s.\n\n',batch.baseline_dv_km_s,r.total_dv_km_s);
fprintf(fid,'35-target witnesses: %d. Maximum witness distance: %.9f km. Minimum altitude: %.9f km.\n\n',r.unique_visit_count,max(r.event_distances_km),r.min_altitude_km);
fprintf(fid,'Two-body only. No J2/ATK validation. Event witnesses certify visits; first-entry times are not exhaustively located.\n\n');
fprintf(fid,'| job | seed | search s | km/s | verified | shots | S2 screened | S2 accepted | burns |\n|---|---|---|---|---|---|---|---|---|\n');
for j=1:n
 fprintf(fid,'| %d | %d | %.2f | %.12f | %d | %d | %d | %d | %d |\n',rows(j,[1 2 3 4 5 6 8 10 13]));
end
end
