function check_v3_connection_reuse()
%CHECK_V3_CONNECTION_REUSE Fresh short controls, never a historical trajectory.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs','v3','preprocessing', ...
 'targets_20260923_release','target_ephemeris.mat'));
c=ctocscreen.v3Defaults(); m=eph.model;
s=struct('schema_version','free_maneuver_v3','dynamics_id','central_j2', ...
 'initial_q',[m.re+600 .0002 .0003 .1 .2 .3], ...
 'maneuver_times_s',[0;500;1200],'duration_s',1800, ...
 'delta_v_km_s',[0 .01 0;0 -.005 .001;.001 -.002 0], ...
 'witness_times_s',nan(35,1),'validation_level','synthetic_test_only');
x=ctocscreen.initialState(s.initial_q,m.mu,m.re).'; anchors=zeros(3,3);
ends=[500;1200;1800];
for k=1:3
 x(4:6)=x(4:6)+s.delta_v_km_s(k,:).';
 x=ctocscreen.v3Arc(x,s.maneuver_times_s(k),ends(k),m,c,false,true);
 anchors(k,:)=x(1:3).';
end
trial=s; trial.delta_v_km_s=trial.delta_v_km_s*1.1;
[result,report]=ctocscreen.v3ConnectControls(trial,anchors,m,c,10);
assert(report.passed&&max(report.endpoint_errors_km)<=c.connection_tolerance_km);
assert(max(abs(result.delta_v_km_s-s.delta_v_km_s),[],'all')<1e-7);
assert(report.reused_integrations>=3);
roots=ctocscreen.v3Roots(eph,ctocscreen.v3Defaults(struct('root_wait_max_s',600)),RandStream('mt19937ar','Seed',888));
assert(roots{1}.seed_wait_s==0&&roots{2}.seed_wait_s>0);
assert(rad2deg(roots{2}.schedule.initial_q(4))>=15);
direct=ctocscreen.v3Roots(eph,ctocscreen.v3Defaults(struct('root_wait_max_s',0)),RandStream('mt19937ar','Seed',888));
for k=1:2
 n=direct{k}; goal=ctocscreen.v3QueryTargets(eph,n.seed_target,n.seed_duration_s);
 h=cross(n.state(1:3),n.state(4:6));
 assert(abs(dot(h,goal))/(norm(h)*norm(goal))<1e-12);
end
root=roots{2}; coast=ctocscreen.v3ApplyAction(root,0,zeros(3,1),root.seed_wait_s,0,eph,c);
assert(coast.t==root.seed_wait_s&&isempty(coast.schedule.maneuver_times_s)&&coast.J==0);
actual=ctocscreen.v3Arc(root.state,0,coast.t,m,c,false,true);
assert(norm(actual(1:3)-coast.state(1:3).')<c.prefix_consistency_km);
[stable,~,checked]=ctocscreen.v3StablePrefix(s,eph,c);
assert(checked.passed&&numel(checked.policy_actions)==3);
freshKeys=ctocscreen.v3ScheduleExperience(stable,eph,c,'prefix_witnesses');
assert(isequal(checked.experience_keys,freshKeys));
smallCap=c; smallCap.search_max_dv_km_s=.01;
cuts=ctocscreen.v3CostReplanCuts(s,smallCap);
assert(any(cuts(1:min(2,numel(cuts)))==0));
table=containers.Map('KeyType','char','ValueType','double');
for k=1:numel(checked.policy_actions)
 a=checked.policy_actions{k};
 feedback=ctocscreen.v3PolicyFeedback(table,a.origin,a.target,a.dt,a.dv,smallCap, ...
  sum(vecnorm(s.delta_v_km_s,2,2)),false);
 assert(feedback.rejected_prefix&&feedback.new_penalty>0);
end
fprintf('CONNECTION_CHECK: fresh three-arc controls recovered; reused integrations=%d.\n',report.reused_integrations);
end
