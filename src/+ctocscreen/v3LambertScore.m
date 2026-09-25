function [score,departure,raw,plane]=v3LambertScore(node,id,dt,eph,c)
%V3LAMBERTSCORE Low-cost seed ranking; no arrival velocity-matching constraint.
goal=ctocscreen.v3QueryTargets(eph,id,node.t+dt);
policy=struct('max_revolutions',ctocscreen.v3LambertRevolutions(node.state(1:3),goal,dt,eph.model.mu,c.guided_max_revolutions),'endpoint_tol_km',.001);
branches=ctocscreen.v3LambertBranches(node.state(1:3),goal,dt,eph.model.mu,policy);
score=Inf; departure=[]; raw=Inf; plane=[]; remaining=node.visited; remaining(id)=true;
[peri,apo]=ctocscreen.v3OsculatingApses(eph.states0(~remaining,:),eph.model.mu);
ranges=[min(apo),max(peri)];
for b=1:numel(branches)
 dv=branches(b).v_depart-node.state(4:6); change=ctocscreen.v3PlaneChange(node.state,dv,c);
 arrival=branches(b).v_arrive;
 [radial,shape]=ctocscreen.v3ContinuationEstimate([goal arrival],remaining,eph,ranges);
 value=norm(dv)+c.plane_penalty_km_s*change.penalty+c.arrival_shape_weight*shape+c.radial_tour_weight*radial;
 if value<score
  score=value; raw=norm(dv); departure=branches(b).v_depart; plane=change;
 end
end
slot=ctocscreen.v3TimeAllowance(node,id,eph);
% Reserve time softly across the mission; do not impose an equal-share leg cap.
score=score+c.cost_time_weight_km_s*dt/max(slot,1);
end
