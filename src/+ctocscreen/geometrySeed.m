function [base,id]=geometrySeed(base,p,stream)
%GEOMETRYSEED Biased starting guess, not a restriction on initial orbit.
id=randi(stream,35); mu=p.mu_km3_s2; a=base.initial_q(1);
h=cross(p.states0(id,1:3),p.states0(id,4:6)); h=h/norm(h);
inc=acos(h(3)); node=atan2(h(1),-h(2));
tof=12000;
for k=1:5
 target=ctocscreen.targetStates(p,id,base.wait_s+tof);
 tof=pi*sqrt(((a+norm(target(1:3)))/2)^3/mu);
end
P=[cos(node) sin(node) 0]; Q=cross(h,P);
direction=-target(1:3)/norm(target(1:3));
u=atan2(dot(direction,Q),dot(direction,P));
u=u-sqrt(mu/a^3)*base.wait_s+.05*randn(stream);
base.initial_q(4:6)=[inc mod(node,2*pi) mod(u,2*pi)];
end
