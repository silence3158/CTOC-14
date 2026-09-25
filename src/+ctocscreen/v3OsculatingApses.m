function [peri,apo]=v3OsculatingApses(states,mu)
%V3OSCULATINGAPSES Instantaneous two-body geometry, not a J2 altitude proof.
radius=vecnorm(states(:,1:3),2,2);
h=cross(states(:,1:3),states(:,4:6),2);
ev=cross(states(:,4:6),h,2)/mu-states(:,1:3)./radius;
eccentricity=vecnorm(ev,2,2);
semimajor=1./(2./radius-sum(states(:,4:6).^2,2)/mu);
peri=semimajor.*(1-eccentricity); apo=semimajor.*(1+eccentricity);
end
