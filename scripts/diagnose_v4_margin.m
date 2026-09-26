function diagnosis=diagnose_v4_margin(pulses,budget)
%DIAGNOSE_V4_MARGIN Why full-history fuel calls stall; run-04 controls are diagnostic input only.
% Hypothesis: after fuel steps, replay fails because active targets or the
% initial orbit land just outside zero-margin acceptance, not because of large misses.
if nargin<1, pulses=12; end
if nargin<2, budget=8; end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
s=load(fullfile(sim,'runs/v4/search/v4_cold_first_600_seed888_04/result.mat'),'result');
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
% Current defaults (run 04 used defaults except budget/seed, which joint ignores).
c=ctocscreen.v4.defaults(); if ~isfield(c,'model_radius_km'), c=s.result.manifest.config; end
node=s.result.best; q=node.q;
q.T=q.tau(pulses+1); q.tau=q.tau(1:pulses); q.u=q.u(1:pulses,:); q.witness(q.witness>q.T)=NaN;
[node.actual,node.trace]=ctocscreen.v4.replay(q,eph,c); node.q=node.trace.q; node.q.witness=node.actual.witness_times_s;
ids=find(node.actual.distance_km<=1); theta=node.actual.witness_times_s(ids);
t0=tic; [child,jr,warm]=ctocscreen.v4.joint(node,ids,theta,'full',node.q.T,eph,c,budget);
elapsed=toc(t0);
last=struct('distance_km',[],'eccentricity',NaN,'height_passed',NaN,'status','');
if ~isempty(warm)
 p=ctocscreen.v4.problem(warm.candidate,ids,theta,'full',node.q.T,eph,c);
 if numel(warm.z)==p.n
  lq=ctocscreen.v4.decode(p,warm.z); [la,~]=ctocscreen.v4.replay(lq,eph,c);
  last=struct('distance_km',la.distance_km(ids),'eccentricity',la.initial_orbit.e, ...
   'sma_km',la.initial_orbit.a,'height_passed',la.height_passed,'status',la.status);
 end
end
childJ=NaN; if ~isempty(child), childJ=child.actual.total_dv_km_s; end
fprintf('prefix pulses=%d active=%d J0=%.9f status=%s it=%d acc=%d rebase=%d final_J=%.9f dJ=%.6f elapsed=%.2f\n', ...
 pulses,numel(ids),node.actual.total_dv_km_s,jr.status,jr.iterations,jr.accepted_steps,jr.rebase_count,childJ, ...
 node.actual.total_dv_km_s-childJ,elapsed);
if ~isempty(last.distance_km)
 fprintf('last iterate: status=%s max active distance=%.6f km, over 0.99: %d, over 1: %d, e=%.9f, height=%d\n', ...
  last.status,max(last.distance_km),sum(last.distance_km>.99),sum(last.distance_km>1),last.eccentricity,last.height_passed);
end
diagnosis=struct('diagnostic_only',true,'pulses',pulses,'ids',ids,'joint',jr,'last',last, ...
 'child_J',childJ,'signature',ctocscreen.v4.signature());
end
