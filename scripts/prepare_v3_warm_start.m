function folder=prepare_v3_warm_start(source,label,budget_s)
%PREPARE_V3_WARM_START Re-aim a historical order/timing in nominal J2.
% Fixed timing is an initializer only. No target velocity/position resets.
if nargin<1, source='runs/screening/beam_v1_s1b_job0001/elite.mat'; end
if nargin<2, label='warm_repair_01'; end
if nargin<3, budget_s=120; end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
old=load(fullfile(sim,source));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
assert(isequal(size(old.p.states0),size(eph.states0))&&max(abs(old.p.states0(:)-eph.states0(:)))<1e-10, ...
 'ctocscreen:v3:warmData','Historical initial target states differ.');
folder=fullfile(sim,'runs/v3/search',label); assert(~isfolder(folder)); mkdir(folder);
c=ctocscreen.v3Defaults(struct('budget_s',budget_s,'shooting_reltol',1e-13,'max_step_s',60));
s=struct('schema_version','free_maneuver_v3','dynamics_id','central_j2', ...
 'initial_q',old.elite.candidate.initial_q,'maneuver_times_s',zeros(0,1), ...
 'delta_v_km_s',zeros(0,3),'duration_s',1,'witness_times_s',nan(35,1),'validation_level','proposal');
x=ctocscreen.initialState(s.initial_q,eph.model.mu,eph.model.re).'; t=0; clock=tic; failure='';
for j=1:numel(old.elite.candidate.order)
 if toc(clock)>=budget_s, failure='budget'; break; end
 try
  depart=old.elite.evaluation.depart_times_s(j); arrive=old.elite.evaluation.arrive_times_s(j);
  id=old.elite.candidate.order(j);
  if depart>t, x=ctocscreen.v3Arc(x,t,depart,eph.model,c,false,true); end
  goal=ctocscreen.v3QueryTargets(eph,id,arrive).';
  dv=ctocscreen.v3GuidedTransfer(x,depart,arrive,goal,eph.model,c,clock);
  [next,~,sol]=ctocscreen.v3Arc([x(1:3);x(4:6)+dv],depart,arrive,eph.model,c,false,true);
  clearance=ctocscreen.v3Height(sol,eph.model,c);
  assert(clearance.passed&&norm(next(1:3)-goal)<=1,'ctocscreen:v3:warmMiss','Rejected corrected leg.');
  s.maneuver_times_s(end+1,1)=depart; s.delta_v_km_s(end+1,:)=dv.';
  s.witness_times_s(id)=arrive; s.duration_s=arrive; x=next; t=arrive;
  fprintf('J2 initializer %d/35: %.9f km/s\n',j,sum(vecnorm(s.delta_v_km_s,2,2)));
 catch err
  failure=[err.identifier ': ' err.message]; break
 end
end
s=ctocscreen.v3Normalize(s,eph.model); verification=ctocscreen.v3Verify(s,eph,c);
if verification.passed, s.validation_level='nominal_j2_independent_verified'; end
elapsed_s=toc(clock); initializer_scope='historical fixed order and timing; J2 re-aimed, free variables released in V3 search';
save(fullfile(folder,'warm_start.mat'),'s','verification','failure','elapsed_s','source','initializer_scope','c');
fprintf('Initializer independent: %d/35 height %d J %.9f time %.3f; %s\n', ...
 verification.visit_count,verification.height_passed,verification.total_dv_km_s,elapsed_s,failure);
end
