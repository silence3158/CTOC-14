function d=v3FeatureDistance(a,b,c)
%V3FEATUREDISTANCE Dimensionless search separation; below one is a near copy.
if strcmp(a.physical_key,b.physical_key), d=0; return; end
s=c.diversity_scales;
d=max([norm(a.initial_state(1:3)-b.initial_state(1:3))/s(1), ...
 norm(a.initial_state(4:6)-b.initial_state(4:6))/s(2), ...
 sqrt(mean(sum((a.states(:,1:3)-b.states(:,1:3)).^2,2)))/s(3), ...
 sqrt(mean(sum((a.states(:,4:6)-b.states(:,4:6)).^2,2)))/s(4), ...
 abs(a.duration_s-b.duration_s)/s(5)]);
if numel(a.maneuver_times_s)~=numel(b.maneuver_times_s)
 d=max(d,2); return;
end
if ~isempty(a.maneuver_times_s)
 d=max([d;abs(a.maneuver_times_s-b.maneuver_times_s)/s(5); ...
  vecnorm(a.delta_v_km_s-b.delta_v_km_s,2,2)/s(6)]);
end
% Only confidently ordered events contribute; near-coincident labels do not.
ta=[a.witness_times_s;a.maneuver_times_s]; tb=[b.witness_times_s;b.maneuver_times_s];
da=ta-ta.'; db=tb-tb.'; tol=c.diversity_event_tolerance_s;
changed=isfinite(da)&isfinite(db)&abs(da)>tol&abs(db)>tol&sign(da)~=sign(db);
if any(changed(:))&&d>0.01, d=max(d,2); end
end
