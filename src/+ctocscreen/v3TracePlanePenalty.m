function penalty=v3TracePlanePenalty(s,trace,c)
%V3TRACEPLANEPENALTY Reuse propagated states, including a terminal impulse.
penalty=0;
for j=1:numel(s.maneuver_times_s)
 index=find(trace.times==s.maneuver_times_s(j),1);
 assert(~isempty(index),'ctocscreen:v3:planeTrace','Missing actual maneuver state.');
 x=trace.states(index,:); x(4:6)=x(4:6)-s.delta_v_km_s(j,:);
 p=ctocscreen.v3PlaneChange(x,s.delta_v_km_s(j,:),c);
 penalty=penalty+p.penalty;
end
end
