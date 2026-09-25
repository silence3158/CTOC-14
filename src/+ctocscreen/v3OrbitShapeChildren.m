function [children,failures]=v3OrbitShapeChildren(node,eph,c)
%V3ORBITSHAPECHILDREN Explicit in-plane coasting preparation, with real impulses.
children={}; failures={}; x=node.state(:); r=norm(x(1:3));
if isfield(node,'transition_kind')&&strcmp(node.transition_kind,'in_plane_orbit_shape')&&node.last_gain==0
 % Full/half corrections already branch; avoid repeated 60-second halvings.
 return;
end
if r<eph.model.re+5000||node.t+60>eph.model.horizon_s, return; end
radial=x(1:3)/r; vr=dot(radial,x(4:6));
if abs(vr)>.8, return; end
transverse=x(4:6)-vr*radial; transverse=transverse/norm(transverse);
dv=sqrt(eph.model.mu/r)*transverse-x(4:6);
if norm(dv)<.03, return; end
for scale=[1 .5]
 impulse=scale*dv;
 if node.J+norm(impulse)>c.search_max_dv_km_s*c.rejected_proposal_factor, continue; end
 try
  child=ctocscreen.v3ApplyAction(node,node.t,impulse,node.t+60,0,eph,c);
  child.transition_kind='in_plane_orbit_shape'; children{end+1}=child;
 catch err
  failures{end+1}=struct('id',err.identifier,'message',err.message);
 end
end
end
