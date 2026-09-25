function [y,P,sol]=v3Arc(x,t0,t1,m,c,variational,independent)
%V3ARC Continuous J2 propagation; no target reset. Absolute times in seconds.
if nargin<6, variational=false; end
if nargin<7, independent=false; end
assert(t0>=0 && t1>=t0 && t1<=m.horizon_s,'ctocscreen:v3:arcTime', ...
 'Invalid arc interval [%.17g, %.17g], horizon %.17g.',t0,t1,m.horizon_s);
x=x(:); P=eye(6); sol=[];
if t1==t0, y=x; return; end
rt=c.shooting_reltol; positionAtol=min(1e-9,10*rt);
at=[positionAtol*ones(3,1);positionAtol*1e-3*ones(3,1)];
if independent, rt=c.verify_reltol; at=[1e-12*ones(3,1);1e-15*ones(3,1)]; end
if variational, x=[x;reshape(P,[],1)]; at=[at;min(1e-11,10*rt)*ones(36,1)]; end
opt=odeset('RelTol',rt,'AbsTol',at,'MaxStep',c.max_step_s);
if independent, sol=ode89(@rhs,[t0 t1],x,opt); else, sol=ode113(@rhs,[t0 t1],x,opt); end
assert(sol.x(end)>=t1,'ctocscreen:v3:integration','Integrator stopped early.');
y=deval(sol,t1); if variational, P=reshape(y(7:end),6,6); end
y=y(1:6); assert(all(isfinite(y)),'ctocscreen:v3:integration','Nonfinite propagation.');
 function d=rhs(t,z)
  if independent
   assert(norm(z(1:3))>m.re/10,'ctocscreen:v3:singularity','Independent trajectory near gravity singularity.');
   a=ctocscreen.v3ReferenceForce(t,z(1:3).',m).'; d=[z(4:6);a];
  elseif variational
   [a,A]=ctocscreen.v3Force(t,z(1:3),m);
   d=[z(4:6);a;reshape([zeros(3),eye(3);A,zeros(3)]*reshape(z(7:end),6,6),[],1)];
  else
   d=[z(4:6);ctocscreen.v3Force(t,z(1:3),m)];
  end
 end
end
