function check_v3_prefix_projection()
%CHECK_V3_PREFIX_PROJECTION Synthetic gate check, never a competition result.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'tests'));
[eph,c,s]=v3SyntheticFixture('separated_visits');
c.shooting_reltol=1e-13; c.joint_seconds=8; c.joint_target_ids=1:8;
s.duration_s=200; s.witness_times_s(9:end)=NaN;
s.maneuver_times_s=100; s.delta_v_km_s=[.001 0 0];
[candidate,info]=ctocscreen.v3JointOptimize(s,eph,c);
assert(strcmp(info.method,'joint_sqp_physical_retraction'));
assert(strcmp(info.status,'prefix_restored')&&~info.replay.passed);
assert(sum(isfinite(candidate.witness_times_s))<35);
assert(ctocscreen.v3RecoveryMerit(info.replay,1:8,1)==0);
assert(info.accepted_steps>=1&&sum(vecnorm(candidate.delta_v_km_s,2,2))<.001);
fprintf('PREFIX_PROJECTION_CHECK: accepted steps=%d, still a partial synthetic trajectory.\n',info.accepted_steps);
end
