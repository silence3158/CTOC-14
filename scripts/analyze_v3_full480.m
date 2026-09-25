function analyze_v3_full480
% Diagnostic replay of the existing run; no search or archive modification.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
a=load(fullfile(sim,'runs/v3/search/v3_full_search_480_01/checkpoint.mat')); state=a.state;
w=load(fullfile(sim,'runs/v3/search/warm_repair_01/warm_start.mat'));
c=state.config; eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
fprintf('Same physical elite: %d\n',strcmp(ctocscreen.v3ScheduleKey(w.s),ctocscreen.v3ScheduleKey(state.elite.schedule)));
fprintf('Joint seconds %.6f; construction seconds %.6f\n',sum(cellfun(@(h)h.diagnostic.elapsed_s,state.history)),sum(cellfun(@(h)h.elapsed_s,state.construction_history)));
p=ctocscreen.v3ShootingProblem(w.s,eph,c); c.joint_iterations=16;
[~,raw,fuel]=ctocscreen.v3SparseFuel(p,p.z0,c,40);
fprintf('Diagnostic fuel: raw %.12f violation %.6g iterations %d\n',fuel.raw_objective,fuel.raw_violation,fuel.output.iterations);
checks=cell(1,5); alphas=[1 .3 .1 .01 .001];
for k=1:numel(alphas)
 s=p.decode(p.z0+alphas(k)*(raw-p.z0)); checks{k}=ctocscreen.v3Replay(s,eph,c,true);
 r=checks{k}; fprintf('alpha %.3g J %.12f visits %d height %d max_distance %.9g\n',alphas(k),r.total_dv_km_s,r.visit_count,r.height_passed,max(r.distance_km));
end
save(fullfile(sim,'runs/v3/development/full480_physical_diagnostic.mat'),'fuel','raw','checks','alphas','c');
end
