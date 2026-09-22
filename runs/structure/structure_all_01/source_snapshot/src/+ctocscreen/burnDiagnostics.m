function tab=burnDiagnostics(source,p)
%BURNDIAGNOSTICS Actual pre/post speeds. Components add in SQUARE, not linearly.
s=source.schedule;r=source.independent;x=r.preburn_states;
before=x(:,4:6);after=before+s.delta_v_km_s;
vb=vecnorm(before,2,2);va=vecnorm(after,2,2);
cosine=sum(before.*after,2)./(vb.*va);cosine=max(-1,min(1,cosine));
direction=sqrt(max(0,2*vb.*va.*(1-cosine)));
tab=table((1:numel(vb))',s.maneuver_times_s,vecnorm(x(:,1:3),2,2)-p.re_km, ...
 vb,va,acos(cosine)*180/pi,va-vb,direction,vecnorm(s.delta_v_km_s,2,2), ...
 'VariableNames',{'burn','time_s','altitude_km','speed_before_km_s','speed_after_km_s', ...
 'velocity_turn_deg','speed_change_km_s','direction_component_km_s','dv_km_s'});
end
