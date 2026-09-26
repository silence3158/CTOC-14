function [value,V,J]=modelMerit(p,z,e,d,lambda,fuel)
%MODELMERIT Eliminate epigraph/slacks for a consistent predicted reduction.
vv=e.visit(:)+e.Jvisit*d;
V=sum(abs(e.eq+e.Jeq*d))+sum(max(0,e.g+e.Jg*d)) ...
 +sum(max(0,vecnorm(reshape(vv,3,p.K),2,1)-e.radius)) ...
 +max(0,norm(e.ecc+e.Jecc*d)-e.ecc_radius);
J=p.sv*sum(vecnorm(reshape(z(p.iu(:))+d(p.iu(:)),3,p.M),2,1));
value=V; if fuel, value=J/p.sv+lambda*V; end
if ~fuel, value=value+p.config.restoration_step_weight*norm(d); end
end
