function report=verify_fragment_final(label)
%VERIFY_FRAGMENT_FINAL Tighter independent replay, targets integrated separately.
root=fileparts(fileparts(mfilename('fullpath'))); cd(root); addpath('src');
folder=fullfile('runs/fragments',label); a=load(fullfile(folder,'elite.mat'));
s=a.elite.schedule; p=a.p; mu=p.mu_km3_s2;
% MATLAB enforces a 100*eps floor; stay above it and record the actual request.
opt=odeset('RelTol',2.3e-14,'AbsTol',1e-15,'MaxStep',600);
x=ctocscreen.initialState(s.initial_q,mu); last=0;
knots=unique([0;s.maneuver_times_s;s.event_times_s;s.duration_s]);
dist=inf(numel(s.event_target_ids),1); minAlt=Inf;
for k=1:numel(knots)
 dt=knots(k)-last; arc=ctocscreen.checkArc(x,dt,p); assert(strcmp(arc.status,'ok'));
 minAlt=min(minAlt,arc.min_altitude_km); x=advance(x,dt);
 ids=find(s.event_times_s==knots(k));
 for j=ids'
  target=advance(p.states0(s.event_target_ids(j),:),knots(k));
  dist(j)=norm(x(1:3)-target(1:3));
 end
 burn=find(s.maneuver_times_s==knots(k));
 if ~isempty(burn), x(4:6)=x(4:6)+s.delta_v_km_s(burn,:); end
 last=knots(k);
end
report=struct('passed',numel(unique(s.event_target_ids(dist<=1)))==35&&minAlt>=200, ...
 'event_distances_km',dist,'max_distance_km',max(dist),'min_altitude_km',minAlt, ...
 'rel_tol',2.3e-14,'abs_tol',1e-15,'max_step_s',600, ...
 'max_distance_change_km',max(abs(dist-a.elite.independent.event_distances_km)), ...
 'method','ode113 fixed impulses; separate target integrations from epoch');
save(fullfile(folder,'tight_verification.mat'),'report'); disp(report);
assert(report.passed);
 function out=advance(in,dt)
  if dt==0, out=in; return; end
  [tt,yy]=ode113(@(~,y)[y(4:6);-mu*y(1:3)/norm(y(1:3))^3],[0 dt],in,opt);
  assert(tt(end)==dt&&all(isfinite(yy(end,:)))); out=yy(end,:);
 end
end
