function [state1,info] = propagateTwoBody(state0,dtSeconds,mu,tol)
%PROPAGATETWOBODY Safeguarded universal Kepler propagation (km, s).
if nargin<3, mu=398600.4415; end
if nargin<4, tol=5e-15; end
validateattributes(state0,{'double'},{'real','finite','numel',6});
validateattributes(dtSeconds,{'double'},{'real','finite','scalar'});
validateattributes(mu,{'double'},{'real','finite','positive','scalar'});
state0=state0(:).'; state1=nan(1,6);
info=struct('status','solver_failure','chi',NaN,'z',NaN,'iterations',0, ...
 'f',NaN,'g',NaN,'fdot',NaN,'gdot',NaN);
r0=state0(1:3); v0=state0(4:6); R=norm(r0);
if R==0, return; end
alpha=2/R-dot(v0,v0)/mu; sm=sqrt(mu); rv=dot(r0,v0)/sm;
dt=dtSeconds;
if alpha>0
 period=2*pi/(sm*alpha^1.5); dt=rem(dt,period);
 if dt>period/2, dt=dt-period; end
 if dt<-period/2, dt=dt+period; end
end
if dt==0
 state1=state0; info.status='ok'; info.chi=0; info.z=0;
 info.f=1; info.g=0; info.fdot=0; info.gdot=1; return;
end
sgn=sign(dt); lo=0; hi=max(1,sqrt(R));
for j=1:100
 [fh,~]=equation(sgn*hi);
 if ~isfinite(fh)||sgn*fh>=0, break; end
 hi=2*hi;
end
x=(lo+hi)/2;
for j=1:160
 [F,dF]=equation(sgn*x); F=sgn*F; info.iterations=j;
 if ~isfinite(F), hi=x; x=(lo+hi)/2; continue; end
 if F>0, hi=x; else, lo=x; end
 if abs(F)<=tol*max(sm*abs(dt),R) || hi-lo<=8*eps(max(1,x))
  info.status='ok'; break;
 end
 trial=x-F/dF;
 if ~isfinite(trial)||trial<=lo||trial>=hi, trial=(lo+hi)/2; end
 x=trial;
end
if ~strcmp(info.status,'ok'), return; end
chi=sgn*x; z=alpha*chi^2; [C,S]=ctocscreen.stumpff(z);
f=1-chi^2*C/R; g=dt-chi^3*S/sm; r=f*r0+g*v0; rmag=norm(r);
if ~isfinite(rmag)||rmag<=0, info.status='solver_failure'; return; end
fdot=sm*(z*S-1)*chi/(rmag*R); gdot=1-chi^2*C/rmag;
state1=[r,fdot*r0+gdot*v0];
if any(~isfinite(state1)), info.status='solver_failure'; state1(:)=NaN; end
info.chi=chi; info.z=z; info.f=f; info.g=g; info.fdot=fdot; info.gdot=gdot;
 function [F,dF]=equation(chi0)
  zz=alpha*chi0^2; [c,s]=ctocscreen.stumpff(zz);
  F=rv*chi0^2*c+(1-R*alpha)*chi0^3*s+R*chi0-sm*dt;
  dF=rv*chi0*(1-zz*s)+(1-R*alpha)*chi0^2*c+R;
 end
end
