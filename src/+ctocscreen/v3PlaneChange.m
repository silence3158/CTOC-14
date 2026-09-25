function detail=v3PlaneChange(state,dv,c)
%V3PLANECHANGE Actual instantaneous osculating inclination and plane rotation.
x=state(:); dv=dv(:); a=cross(x(1:3),x(4:6)); b=cross(x(1:3),x(4:6)+dv);
if norm(a)<eps||norm(b)<eps
 detail=struct('inclination_change_deg',180,'plane_rotation_deg',180,'penalty',Inf); return;
end
a=a/norm(a); b=b/norm(b);
ia=acosd(max(-1,min(1,a(3)))); ib=acosd(max(-1,min(1,b(3))));
change=abs(ib-ia); angle=acosd(max(-1,min(1,dot(a,b))));
detail=struct('inclination_change_deg',change,'plane_rotation_deg',angle, ...
 'penalty',(max(0,change-c.plane_change_threshold_deg)/c.plane_change_threshold_deg)^2);
end
