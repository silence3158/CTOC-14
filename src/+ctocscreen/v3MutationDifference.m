function report=v3MutationDifference(a,b,c,m)
%V3MUTATIONDIFFERENCE Reject numerical return to the parent, not distant families.
tol=c.mutation_difference_tolerances;
xa=ctocscreen.initialState(a.initial_q,m.mu,m.re); xb=ctocscreen.initialState(b.initial_q,m.mu,m.re);
ka=vecnorm(a.delta_v_km_s,2,2)>tol(2); kb=vecnorm(b.delta_v_km_s,2,2)>tol(2);
ta=a.maneuver_times_s(ka); tb=b.maneuver_times_s(kb);
da=a.delta_v_km_s(ka,:); db=b.delta_v_km_s(kb,:);
report=struct('initial_position_km',norm(xa(1:3)-xb(1:3)), ...
 'initial_velocity_km_s',norm(xa(4:6)-xb(4:6)), ...
 'effective_maneuver_count_changed',numel(ta)~=numel(tb), ...
 'time_change_s',abs(a.duration_s-b.duration_s),'pulse_change_km_s',0,'novel',false);
if numel(ta)==numel(tb)&&~isempty(ta)
 report.time_change_s=max(report.time_change_s,max(abs(ta-tb)));
 report.pulse_change_km_s=max(vecnorm(da-db,2,2));
end
report.novel=report.effective_maneuver_count_changed||report.initial_position_km>tol(1) ...
 ||report.initial_velocity_km_s>tol(2)||report.time_change_s>tol(3)||report.pulse_change_km_s>tol(2);
end
