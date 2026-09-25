function check_v3_route_guidance()
%CHECK_V3_ROUTE_GUIDANCE Focused fallback and radial-order sanity checks.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs','v3','preprocessing', ...
 'targets_20260923_release','target_ephemeris.mat'));
c=ctocscreen.v3Defaults(struct('budget_s',10)); m=eph.model;
x=[7000;0;0;0;sqrt(m.mu/7000);0]; dt=1500;
goal=ctocscreen.v3Arc(x,0,dt,m,c);
clock=tic;
[dv,stats,alternatives]=ctocscreen.v3GuidedTransfer(x,0,dt,goal(1:3),m,c,clock,[0;1;0]);
assert(norm(dv)<.001&&stats.branches_tried>1&&~isempty(alternatives));
assert(ctocscreen.v3ContinuationEstimate(x,true(35,1),eph)==0);
visited=true(35,1); visited([1 12])=false;
middle=ctocscreen.v3ContinuationEstimate([26562 0 0 0 sqrt(m.mu/26562) 0],visited,eph);
inside=ctocscreen.v3ContinuationEstimate(x,visited,eph);
assert(inside>middle&&isfinite(middle));
rp=7000; ra=45000; r=26562; h=sqrt(2*m.mu*rp*ra/(rp+ra));
vt=h/r; vr=sqrt(m.mu*(2/r-2/(rp+ra))-vt^2);
[accessible,shape]=ctocscreen.v3ContinuationEstimate([r 0 0 vr vt 0],visited,eph);
assert(accessible<1e-10&&shape==0&&middle>accessible);
fprintf('ROUTE_CHECK: failed supplied branch retries; radial surrogate finite and empty-task zero.\n');
end
