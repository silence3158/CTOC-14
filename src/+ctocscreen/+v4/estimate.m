function [H,detail]=estimate(node,eph,c,stream)
%ESTIMATE Remaining-cost heuristic H (km/s) for beam ranking; not a bound.
% Uses the cheapest crossing-event Lambert proposal per remaining target
% from this node (eval_v4_estimator: Pearson 0.985 with realized legs, the
% two-body tangential estimate 0.261). Targets without a proposal in the
% window take the largest found value. H = weight * sum of the cheapest
% fraction of per-target values.
remaining=find(node.actual.distance_km>1);
detail=struct('events',0,'lambert',0,'seconds',0,'per_target',nan(35,1));
if isempty(remaining), H=0; return; end
clock=tic;
cc=c; cc.event_candidates=max(c.event_candidates,2*numel(remaining)); cc.events_per_target=1; cc.departures_per_event=2;
[ev,info]=ctocscreen.v4.events(node,eph,cc);
[~,cr]=ctocscreen.v4.crossing(node,eph,cc,ev,info,stream,c.estimate_seconds);
v=cr.target_best(remaining); ok=isfinite(v);
if ~any(ok), H=c.no_event_cost_km_s*numel(remaining);
else
 v(~ok)=max(v(ok)); s=sort(v); n=max(1,ceil(c.heuristic_fraction*numel(s)));
 H=c.heuristic_weight*sum(s(1:n))*numel(s)/n;
end
detail.events=numel(ev); detail.lambert=cr.lambert; detail.seconds=toc(clock); detail.per_target(remaining)=v;
end
