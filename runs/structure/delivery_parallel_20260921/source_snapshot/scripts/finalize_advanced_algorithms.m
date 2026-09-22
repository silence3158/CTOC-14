function finalize_advanced_algorithms()
%FINALIZE_ADVANCED_ALGORITHMS Verify nonwinning demonstrations and export best.
root=fileparts(fileparts(mfilename('fullpath')));cd(root);addpath('src');
folder='runs/fragments/algorithms_20260921';a=load(fullfile(folder,'study.mat'));study=a.study;
base=load('runs/fragments/pilot8min_20260921/elite.mat');p=base.p;elite=base.elite;
joint=load(fullfile(folder,'joint_corrected.mat'));old=study.outputs{1};
periodic=load(fullfile(folder,'joint_periodic.mat'));
warmPeriodic=load(fullfile(folder,'joint_periodic_warm.mat'));
for candidate={old,joint.o,periodic.o,warmPeriodic.o}
 o=candidate{1};if o.independent.passed&&o.independent.total_dv_km_s<elite.independent.total_dv_km_s-1e-8,elite=o;end
end
cold=study.outputs{2};vals=cellfun(@(x)x.evaluation.total_dv_km_s,cold.complete);
[~,j]=min(vals);cold_check=ctocscreen.propagateSchedule(cold.complete{j}.schedule,p,true);
multi=study.outputs{3};fragment_checks=cell(1,numel(multi.fragments));
cfg=ctocscreen.advancedDefaults();
for j=1:numel(multi.fragments)
 f=multi.fragments{j};entry=find(abs(base.elite.schedule.event_times_s-f.t_in_s)<1e-8,1);
 if isempty(entry),entry=0;end
 [s,r]=ctocscreen.appendFragmentSuffix(base.elite.schedule,entry,f,p,cfg);
 ir=struct('passed',false);if r.passed,ir=ctocscreen.propagateSchedule(s,p,true);end
 fragment_checks{j}=struct('counts',f.counts,'ids',f.event_target_ids,'screened',r,'independent',ir);
end
save(fullfile(folder,'demonstration_checks.mat'),'cold_check','fragment_checks','-v7.3');
study.final_best_verified_dv=elite.independent.total_dv_km_s;
study.corrected_joint=joint.o;save(fullfile(folder,'study.mat'),'study','-v7.3');
study.periodic_joint=periodic.o;save(fullfile(folder,'study.mat'),'study','-v7.3');
study.warm_periodic_joint=warmPeriodic.o;save(fullfile(folder,'study.mat'),'study','-v7.3');
cfg=elite.config;save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
s=elite.schedule;r=elite.independent;
writetable(table(s.event_target_ids,s.event_times_s,r.event_distances_km, ...
 'VariableNames',{'target','time_s','distance_km'}),fullfile(folder,'events.csv'));
writetable(table(s.maneuver_times_s,s.delta_v_km_s(:,1),s.delta_v_km_s(:,2),s.delta_v_km_s(:,3), ...
 vecnorm(s.delta_v_km_s,2,2),'VariableNames',{'time_s','dvx','dvy','dvz','dv_km_s'}),fullfile(folder,'maneuvers.csv'));
f=figure('Visible','off','Color','w');plot(old.history(:,1),old.history(:,2),'DisplayName','Initial joint implementation');hold on;
plot(joint.o.history(:,1),joint.o.history(:,2),'DisplayName','Smooth derivative evaluation');
plot(periodic.o.history(:,1),periodic.o.history(:,2),'DisplayName','Unwrapped periodic angles');
plot(warmPeriodic.o.history(:,1),warmPeriodic.o.history(:,2),'DisplayName','Unwrapped angles, warm start');
yline(study.baseline,'--','DisplayName','Previous verified best');legend('Location','best');grid on;
xlabel('Optimization time (s)');ylabel('Screened total delta-V (km/s)');
exportgraphics(f,fullfile(folder,'joint_convergence.png'),'Resolution',160);close(f);
fprintf('FINAL_ADVANCED best%.12f gain%.6f%% days%.9f max_m%.9f altitude%.9f\n', ...
 r.total_dv_km_s,100*(1-r.total_dv_km_s/study.baseline),s.duration_s/86400,max(r.event_distances_km)*1000,r.min_altitude_km);
fprintf('COLD_CHECK passed%d dv%.12f\n',cold_check.passed,cold_check.total_dv_km_s);
for j=1:numel(fragment_checks)
 x=fragment_checks{j};fprintf('MULTI_CHECK %s screened%d independent%d dv%.9f\n',mat2str(x.ids'),x.screened.passed,x.independent.passed,x.screened.total_dv_km_s);
end
disp(elite.plan.initial_q);disp(joint.o.output);disp(joint.o.block_output);
end
