function slot=v3TimeAllowance(node,id,eph)
%V3TIMEALLOWANCE Soft time reservation scaled by unvisited orbital periods.
x=eph.states0; r=vecnorm(x(:,1:3),2,2);
inverseA=2./r-sum(x(:,4:6).^2,2)/eph.model.mu;
a=r; bound=inverseA>0; a(bound)=1./inverseA(bound);
period=2*pi*sqrt(a.^3/eph.model.mu);
left=find(~node.visited);
slot=(eph.model.horizon_s-node.t)*period(id)/sum(period(left));
end
