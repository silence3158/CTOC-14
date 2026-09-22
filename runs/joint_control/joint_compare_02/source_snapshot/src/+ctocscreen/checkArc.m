function result = checkArc(state0,dtSeconds,problem,~)
%CHECKARC Analytic radial extrema on a finite forward two-body conic.
result=struct('status','solver_failure','method','conic_extrema', ...
 'min_radius_km',NaN,'min_altitude_km',NaN,'endpoint_norm_km',NaN, ...
 'propagation_status','solver_failure','periapsis_time_s',NaN);
if ~isscalar(dtSeconds)||~isfinite(dtSeconds)||dtSeconds<0, return; end
mu=problem.mu_km3_s2; re=problem.re_km;
[s1,inf1]=ctocscreen.propagateTwoBody(state0,dtSeconds,mu);
result.propagation_status=inf1.status;
if ~strcmp(inf1.status,'ok'), return; end
r=state0(1:3); v=state0(4:6); R=norm(r); h=cross(r,v); H=norm(h);
if H<1e-10, result.method='radial_degenerate'; return; end
evec=cross(v,h)/mu-r/R; ecc=norm(evec); p=H^2/mu;
alpha=2/R-dot(v,v)/mu; tperi=Inf;
if ecc<1e-10
 result.method='circular_analytic';
else
 cosnu=max(-1,min(1,dot(evec,r)/(ecc*R)));
 sinnu=dot(r,v)*H/(mu*ecc*R);
 nu=atan2(sinnu,cosnu);
 if alpha>1e-12
  a=1/alpha;
  E=atan2(sqrt(max(0,1-ecc^2))*sin(nu),ecc+cos(nu));
  M=mod(E-ecc*sin(E),2*pi); n=sqrt(mu/a^3);
  if M<1e-12, tperi=0; else, tperi=(2*pi-M)/n; end
  result.method='elliptic_analytic';
 elseif alpha < -1e-12
  A=-1/alpha; F=asinh(dot(r,v)/(ecc*sqrt(mu*A)));
  tperi=-(ecc*sinh(F)-F)/sqrt(mu/A^3);
  result.method='hyperbolic_analytic';
 else
  rp=p/(1+ecc);
  if abs(alpha)>1e-16
   if alpha>0
    E=atan2(sqrt(max(0,1-ecc^2))*sin(nu),ecc+cos(nu));
    chi=E/sqrt(alpha);
   else
    F=asinh(dot(r,v)*sqrt(-alpha)/(ecc*sqrt(mu)));
    chi=F/sqrt(-alpha);
   end
   [~,S]=ctocscreen.stumpff(alpha*chi^2);
   tperi=-((1-rp*alpha)*chi^3*S+rp*chi)/sqrt(mu);
   if alpha>0 && tperi<0, tperi=tperi+2*pi/sqrt(mu*alpha^3); end
  else
   D=tan(nu/2); tperi=-0.5*sqrt(p^3/mu)*(D+D^3/3);
  end
  result.method='near_parabolic_analytic';
 end
end
minimum=min(R,norm(s1(1:3)));
if tperi>=-1e-8 && tperi<=dtSeconds+1e-8, minimum=min(minimum,p/(1+ecc)); end
result.status='ok'; result.min_radius_km=minimum;
result.min_altitude_km=minimum-re; result.endpoint_norm_km=norm(s1(1:3));
result.periapsis_time_s=tperi;
end
