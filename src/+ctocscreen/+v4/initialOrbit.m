function o=initialOrbit(x,m)
%INITIALORBIT Osculating initial-orbit constraints and exact Jacobians.
r=x(1:3); v=x(4:6); r=r(:); v=v(:); R=norm(r);
o.alpha=2/R-dot(v,v)/m.mu;
o.evec=(dot(v,v)*r-dot(r,v)*v)/m.mu-r/R;
o.a=1/o.alpha; o.e=norm(o.evec);
o.alpha_jac=[-2*r.'/R^3,-2*v.'/m.mu];
Er=(dot(v,v)*eye(3)-v*v.')/m.mu-eye(3)/R+r*r.'/R^3;
Ev=(2*r*v.'-v*r.'-dot(r,v)*eye(3))/m.mu;
o.ecc_jac=[Er,Ev];
o.passed=isfinite(o.a)&&o.a>=m.re+590&&o.a<=m.re+610&&o.e<.001;
end
