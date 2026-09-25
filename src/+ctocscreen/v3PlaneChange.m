function detail=v3PlaneChange(state,dv,c)
%V3PLANECHANGE Actual instantaneous osculating inclination and plane rotation.
x=state(:); dv=dv(:); a=cross(x(1:3),x(4:6)); b=cross(x(1:3),x(4:6)+dv);
if norm(a)<eps||norm(b)<eps
 detail=struct('inclination_change_deg',180,'plane_rotation_deg',180,'penalty',Inf); return;
end
a=a/norm(a); b=b/norm(b);
ia=acosd(max(-1,min(1,a(3)))); ib=acosd(max(-1,min(1,b(3))));
change=abs(ib-ia); angle=acosd(max(-1,min(1,dot(a,b))));
% Monotone and saturating: a 5 deg change scores 0.5, 10 deg scores 1.0, and
% the score approaches the saturation value instead of being clipped, so
% "slightly over threshold" and "catastrophic plane change" stay distinguishable.
excess=max(0,change-c.plane_change_threshold_deg);
capacity=max(eps,c.plane_penalty_saturation-1);
detail=struct('inclination_change_deg',change,'plane_rotation_deg',angle, ...
 'penalty',c.plane_penalty_saturation*excess/(excess+capacity*c.plane_change_threshold_deg));
end
