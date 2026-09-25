function a = v3J2Acceleration(t,r,model)
%V3J2ACCELERATION Central + degree-2 zonal gravity in native assumed GCRF.
% r: N-by-3 km; t: scalar or N-by-1 seconds. Output N-by-3 km/s^2.
validateattributes(r,{'double'},{'2d','ncols',3,'real','finite'});
validateattributes(t,{'double'},{'vector','real','finite','nonempty'});
t=t(:); n=size(r,1);
assert(isscalar(t)||numel(t)==n,'ctocscreen:v3:shape','Time/state sizes disagree.');
assert(all(t>=0 & t<=model.horizon_s),'ctocscreen:v3:timeRange','Time outside model range.');
% Uniform pole grid: constant-time indexing, avoiding interp1 in ODE calls.
h=model.pole_times(2)-model.pole_times(1);
idx=min(floor(t/h)+1,numel(model.pole_times)-1); u=(t-model.pole_times(idx))/h;
k=model.pole(idx,:).*(1-u)+model.pole(idx+1,:).*u;
k=k./sqrt(sum(k.^2,2));
rho2=sum(r.^2,2); assert(all(rho2>0),'ctocscreen:v3:zeroRadius','Zero radius.');
z=sum(r.*k,2); q=z.^2./rho2;
a=-model.mu*r./rho2.^1.5 + 1.5*model.j2*model.mu*model.re^2./rho2.^2.5 .* ...
 ((5*q-1).*r-2*z.*k);
end
