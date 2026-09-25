function count=v3LambertRevolutions(r1,r2,dt,mu,limit)
%V3LAMBERTREVOLUTIONS Necessary two-body bound, cf. Izzo T/pi in lamberthub.
% This only skips impossible Lambert seed intervals; it is not a J2 limit.
s=(norm(r1)+norm(r2)+norm(r2(:)-r1(:)))/2;
T=sqrt(2*mu/s^3)*dt;
count=min(limit,max(0,floor(T/pi+1e-12)));
end
