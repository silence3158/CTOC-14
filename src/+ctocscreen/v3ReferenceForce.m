function a=v3ReferenceForce(t,r,m)
%V3REFERENCEFORCE Separately arranged zonal formula; N-by-3 states.
h=m.pole_times(2); j=min(floor(t/h)+1,size(m.pole,1)-1);
u=(t-m.pole_times(j))/h; k=(1-u)*m.pole(j,:)+u*m.pole(j+1,:); k=k/norm(k);
R=vecnorm(r,2,2); axial=(r*k.')*k; transverse=r-axial;
q=(r*k.').^2./R.^2; f=1.5*m.j2*(m.re./R).^2;
a=-m.mu./R.^3.*((1+f.*(1-5*q)).*transverse+(1+f.*(3-5*q)).*axial);
end
