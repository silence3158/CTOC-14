function out=v3Height(sol,m,c)
%V3HEIGHT Adaptive interval lower bounds, not endpoint-only clearance.
% Bootstrap: before the first R=Rmin crossing acceleration is bounded by A.
% A midpoint cannot reach Rmin if |r_mid|-|v_mid|h-Ah^2/2 exceeds it.
% Numerical integration uncertainty is handled by an explicit margin.
Rmin=m.re+200; A=m.mu/Rmin^2*(1+6*abs(m.j2)*(m.re/Rmin)^2);
stack=[sol.x(1),sol.x(end)]; lower=Inf; sampled=Inf; checked=0; out.passed=true;
while ~isempty(stack)
 ab=stack(end,:); stack(end,:)=[]; mid=mean(ab); h=diff(ab)/2;
 y=deval(sol,[ab(1),mid,ab(2)]); rr=vecnorm(y(1:3,:),2,1);
 sampled=min(sampled,min(rr)); checked=checked+1;
 if min(rr)<Rmin, out.passed=false; break; end
 bound=rr(2)-norm(y(4:6,2))*h-0.5*A*h*h;
 if bound>=Rmin+c.height_margin_km
  lower=min(lower,bound); continue
 end
 if diff(ab)<=c.height_min_interval_s
  out.passed=false; break
 end
 stack=[stack;ab(1),mid;mid,ab(2)]; %#ok<AGROW>
end
out.lower_bound_altitude_km=lower-m.re;
out.sampled_min_altitude_km=sampled-m.re; out.intervals=checked;
out.method='adaptive acceleration-bound intervals with numerical margin';
end
