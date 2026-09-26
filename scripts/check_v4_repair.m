function report=check_v4_repair()
%CHECK_V4_REPAIR Diagnostic controls are never inputs to the cold search runner.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
s=load(fullfile(sim,'runs/v4/search/v4_cold_first_300_seed888_01/result.mat'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(); parent=s.result.best; q=parent.q;
q.T=q.tau(6); q.tau=q.tau(1:5); q.u=q.u(1:5,:); q.witness(q.witness>q.T)=NaN;
[parent.actual,parent.trace]=ctocscreen.v4.replay(q,eph,c); parent.q=parent.trace.q;
ids=find(parent.actual.distance_km<=c.search_radius_km);
[child,jr,warm]=ctocscreen.v4.joint(parent,ids,parent.actual.witness_times_s(ids),'full',q.T,eph,c,12);
assert(~isempty(child)&&jr.accepted_steps>0,'No accepted nonlinear step on real prefix.');
p=ctocscreen.v4.problem(parent,ids,parent.actual.witness_times_s(ids),'full',q.T,eph,c);
lastq=ctocscreen.v4.decode(p,warm.z);
[lastActual,~]=ctocscreen.v4.replay(lastq,eph,c);
diagnostic=struct('joint',jr,'last_actual',lastActual,'warm',rmfield(warm,{'candidate','best'}));
name=char(datetime('now','Format','yyyyMMdd_HHmmss'));
save(fullfile(sim,'runs/v4/development',['repair_attempt_' name '.mat']),'diagnostic');
disp(jr); disp([ids,lastActual.distance_km(ids)]);
fprintf('Internal J %.12f; incumbent J %.12f\n',lastActual.total_dv_km_s,child.actual.total_dv_km_s);
assert(jr.active_passed&&child.actual.total_dv_km_s<parent.actual.total_dv_km_s-1e-7, ...
 'No actual cost improvement retaining all active real targets.');
[independent,~]=ctocscreen.v4.replay(child.q,eph,c,true);
assert(all(independent.distance_km(ids)<=1)&&independent.height_passed&&independent.initial_passed);
[rebuilt,rr]=ctocscreen.v4.rebuild(parent,eph,c,1);
assert(~isempty(rebuilt)&&rebuilt.q.T<parent.q.T&&numel(rebuilt.q.tau)<numel(parent.q.tau));
assert(rebuilt.actual.visit_count<parent.actual.visit_count,'Revoked suffix must lose coverage.');
assert(numel(rebuilt.actual.distance_km)==35,'Rebuilding must retain the entire mission.');
report=struct('diagnostic_only',true,'prefix_targets',ids,'joint',jr,'rebuild',rr, ...
 'before_dv',parent.actual.total_dv_km_s,'after_dv',child.actual.total_dv_km_s, ...
 'independent',independent,'warm_trust',warm.trust,'signature',ctocscreen.v4.signature());
name=char(datetime('now','Format','yyyyMMdd_HHmmss'));
save(fullfile(sim,'runs/v4/development',['repair_checks_' name '.mat']),'report');
fprintf('V4 REPAIR CHECK: actual prefix %d targets, %.12f -> %.12f km/s; accepted=%d; independent=%d\n', ...
 numel(ids),report.before_dv,report.after_dv,jr.accepted_steps,sum(independent.distance_km(ids)<=1));
fprintf('V4 REBUILD CHECK: %d -> %d visits, retained full 35-target mission\n',parent.actual.visit_count,rebuilt.actual.visit_count);
end
