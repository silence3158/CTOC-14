function report=check_v4_crossing_t1(label)
%CHECK_V4_CROSSING_T1 Do crossing-event proposals give cheaper next legs?
% Compares two-body Lambert proposal Delta-V from the SAME nodes: cold roots
% (seed 888) and one mid-mission node from run-04 controls (diagnostic input
% only, never a search input). Distinguishes: H1 crossing proposals lower the
% cheapest and median next-leg cost vs the old any-time grid; H0 no change.
% Proposal costs are two-body seeds, not J2-corrected or mission results.
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(); stream=RandStream('mt19937ar','Seed',888);
nodes={}; for k=1:4, nodes{end+1}=ctocscreen.v4.root(k,eph,c,stream); end %#ok<AGROW>
s=load(fullfile(sim,'runs/v4/search/v4_cold_first_600_seed888_04/result.mat'),'result');
for cut=[6 12]
 q=s.result.best.q; q.T=q.tau(cut+1); q.tau=q.tau(1:cut); q.u=q.u(1:cut,:); q.witness(q.witness>q.T)=NaN;
 [a,tr]=ctocscreen.v4.replay(q,eph,c); n=nodes{1}; n.q=tr.q; n.q.witness=a.witness_times_s; n.actual=a; n.trace=tr;
 nodes{end+1}=n; %#ok<AGROW>
end
rows=zeros(numel(nodes),8);
for k=1:numel(nodes)
 n=nodes{k};
 t0=tic; [ev,info]=ctocscreen.v4.events(n,eph,c); [cand,cr]=ctocscreen.v4.crossing(n,eph,c,ev,info,stream,30); tc=toc(t0);
 t0=tic; old=oldGrid(n,eph,c,stream); to=toc(t0);
 a=sort([cand.dv]); b=sort(old);
 rows(k,:)=[n.actual.visit_count,numel(a),min(a),median(a(1:min(end,20))),numel(b),min(b),median(b(1:min(end,20))),tc-to];
 fprintf('node %d visits=%2d | crossing n=%3d min=%.3f med20=%.3f (%.2fs, %d Lambert) | old n=%3d min=%.3f med20=%.3f (%.2fs)\n', ...
  k,rows(k,1),rows(k,2),rows(k,3),rows(k,4),tc,cr.lambert,rows(k,5),rows(k,6),rows(k,7),to);
end
report=struct('diagnostic_only',true,'rows',rows,'columns',{{'visits','cross_n','cross_min','cross_med20','old_n','old_min','old_med20','time_diff'}}, ...
 'signature',ctocscreen.v4.signature());
save(fullfile(sim,'runs/v4/development',['crossing_t1_' label '.mat']),'report');
end
function dv=oldGrid(parent,eph,c,stream)
% The pre-2026-09-26 expand.m target/time grid, all pairs, no time cut.
m=eph.model; q=parent.q; t0=q.T; x=parent.actual.final_state; dv=[];
remaining=find(parent.actual.distance_km>1); available=m.horizon_s-t0;
look=min([available,c.lookahead_s,3*available/max(1,numel(remaining))]);
dt=unique([max(600,look*[.06 .12 .25 .45 .7 1]),min(look,parent.seed_duration)]);
dt=dt(dt>0&dt<=available); if numel(dt)>c.time_samples, dt=dt(round(linspace(1,numel(dt),c.time_samples))); end
h=cross(x(1:3),x(4:6)); h=h/max(norm(h),eps);
rr=ctocscreen.v3QueryTargets(eph,remaining,t0+min(look,median(dt)));
[~,rank]=sort(abs(rr*h)+.05*abs(vecnorm(rr,2,2)-norm(x(1:3))));
chosen=remaining(rank(1:min(c.target_count,numel(rank))));
policy=struct('max_revolutions',c.max_revolutions,'endpoint_tol_km',.005);
for id=chosen(:).'
 for d=dt
  rt=ctocscreen.v3QueryTargets(eph,id,t0+d);
  try, b=ctocscreen.v3LambertBranches(x(1:3),rt,d,m.mu,policy); catch, continue; end
  for j=1:numel(b), dv(end+1)=norm(b(j).v_depart(:)-x(4:6)); end %#ok<AGROW>
 end
end
end
