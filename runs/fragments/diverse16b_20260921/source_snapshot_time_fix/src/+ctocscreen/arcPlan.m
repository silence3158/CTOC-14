function plan=arcPlan(s,r)
%ARCPLAN Encode burn arcs with one or more visits; reject zero-visit burns.
M=numel(s.maneuver_times_s); ends=zeros(M,1); starts=zeros(M,1); waits=zeros(M,1);
for m=1:M
 upper=Inf;if m<M,upper=s.maneuver_times_s(m+1);end
 % Visits at a burn epoch are preburn witnesses, hence belong to the preceding arc.
 ix=find(s.event_times_s>s.maneuver_times_s(m)+1e-7&s.event_times_s<=upper+1e-7);
 assert(~isempty(ix),'Zero-visit arcs require the generic direct-shooting operator.');
 starts(m)=ix(1);ends(m)=ix(end);prev=0;if m>1,prev=s.event_times_s(ends(m-1));end
 waits(m)=s.maneuver_times_s(m)-prev;
 if waits(m)<0&&waits(m)>=-1e-7,waits(m)=0;end
end
assert(starts(1)==1&&ends(end)==numel(s.event_times_s));
plan=struct('initial_q',s.initial_q,'ids',s.event_target_ids(:),'times',s.event_times_s(:), ...
 'counts',ends-starts+1,'waits',waits,'reference_v',r.preburn_states(:,4:6)+s.delta_v_km_s);
end
