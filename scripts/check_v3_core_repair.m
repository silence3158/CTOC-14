function check_v3_core_repair
% Development-only complete-task check; not a formal search batch.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
w=load(fullfile(sim,'runs/v3/search/warm_repair_01/warm_start.mat'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v3Defaults(struct('joint_seconds',60,'shooting_reltol',1e-12,'max_step_s',60));
[s,info]=ctocscreen.v3JointOptimize(w.s,eph,c);
save(fullfile(sim,'runs/v3/development/core_repair_check.mat'),'s','info','c');
fprintf('FULL J %.12f visits %d accepted_steps %d elapsed %.3f reason %s\n', ...
 info.replay.total_dv_km_s,info.replay.visit_count,info.accepted_steps,info.elapsed_s,info.failure_reason);
for k=1:numel(info.trials)
 r=info.trials{k}; fprintf('trial %d alpha %.3g radius %.3g connect %d J %.12f accept %d reason %s\n', ...
 k,r.alpha,r.radius,r.connection.passed,r.actual_dv_km_s,r.accepted,r.connection.reason);
 if ~isempty(r.verification), fprintf('visits %d maxdist %.9g\n',r.verification.visit_count,max(r.verification.distance_km)); end
end
end
