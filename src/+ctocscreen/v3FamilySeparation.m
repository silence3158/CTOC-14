function [separation,detail]=v3FamilySeparation(t,dt,eph,id,c,homeBand,homeInc)
%V3FAMILYSEPARATION Ranking-only cost of committing the tour to another family.
% Altitude-band and plane changes cannot be undone for free: a tour that hops
% between the GEO, MEO and inner bands, or between the 0/29/43/55/63 degree
% planes, pays a round trip for each hop. This estimates that commitment so the
% ranking finishes a family before leaving it. It is NOT a physical constraint,
% NOT a lower bound, and never enters the reported delta-V.
% The target state must be taken at the ARRIVAL epoch t+dt: the target moves
% along its orbit, so using the departure epoch overstates the plane/band gap.
[r,v]=ctocscreen.v3QueryTargets(eph,id,t+dt);
r=r(1,:); v=v(1,:);
radius=norm(r);
invA=2/radius-sum(v.^2)/eph.model.mu;
a=1/max(invA,eps);
h=cross(r,v);
inc=acosd(max(-1,min(1,h(3)/max(norm(h),eps))));
bandChange=abs(round(a/2.5e3)-round(homeBand))>0;
bandCost=double(bandChange)*c.band_separation_km_s;
incDiff=abs(inc-homeInc)*c.plane_separation_weight;
separation=bandCost+incDiff;
detail=struct('target_band_km',a,'target_inc_deg',inc,'band_change',bandChange, ...
 'band_cost_km_s',bandCost,'plane_cost_km_s',incDiff,'arrival_epoch_s',t+dt);
end
