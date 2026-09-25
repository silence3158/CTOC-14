function [a,A]=v3Force(t,r,m)
%V3FORCE J2 acceleration and analytic spatial Jacobian for variational ODE.
r=r(:); h=m.pole_times(2); j=min(floor(t/h)+1,size(m.pole,1)-1);
u=(t-m.pole_times(j))/h; k=((1-u)*m.pole(j,:)+u*m.pole(j+1,:)).'; k=k/norm(k);
R=norm(r); assert(R>m.re/10,'ctocscreen:v3:singularity','Trajectory near gravity singularity.');
z=k.'*r; s=5*z*z/R^2-1; b=s*r-2*z*k; C=1.5*m.j2*m.mu*m.re^2;
a=-m.mu*r/R^3+C*b/R^5;
if nargout>1
 ds=10*z/R^2*k.'-10*z*z/R^4*r.';
 A=-m.mu*(eye(3)/R^3-3*(r*r.')/R^5) ...
   +C*((s*eye(3)+r*ds-2*(k*k.'))/R^5-5*(b*r.')/R^7);
end
end
