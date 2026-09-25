function [cuts,retained]=v3CostReplanCuts(s,c)
%V3COSTREPLANCUTS Rebuild before overspending; reserve room for remaining visits.
dv=vecnorm(s.delta_v_km_s,2,2); spent=[0;cumsum(dv(1:end-1))];
cuts=[]; retained=[];
if isempty(dv)||sum(dv)<=c.search_max_dv_km_s, return; end
cap=c.search_max_dv_km_s;
% The first burn is a launch investment. Subsequent budget pacing is only a
% proposal heuristic, never a proof that another prefix cannot be improved.
visits=arrayfun(@(t)sum(s.witness_times_s<=t),s.maneuver_times_s);
allowance=min(.9*cap,dv(1)+(cap-dv(1))*visits/35);
crossing=find(cumsum(dv)>allowance&spent<cap,1);
if ~isempty(crossing), cuts=s.maneuver_times_s(crossing); end
for fraction=[.5 .75 .25 .9]
 j=find(spent<=fraction*cap,1,'last');
 if ~isempty(j), cuts(end+1)=s.maneuver_times_s(j); end
end
if isempty(cuts), cuts=0; else, cuts=unique([cuts(1) 0 cuts(2:end)],'stable'); end
retained=arrayfun(@(t)sum(dv(s.maneuver_times_s<t)),cuts);
keep=retained<cap; cuts=cuts(keep); retained=retained(keep);
end
