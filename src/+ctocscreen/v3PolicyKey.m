function key=v3PolicyKey(node,kind,target,dt)
%V3POLICYKEY Coarse state-conditioned backoff; never replaces exact state keys.
x=node.state(:); r=x(1:3); v=x(4:6); h=cross(r,v); h=h/max(norm(h),eps);
inc=acosd(max(-1,min(1,h(3)))); raan=atan2d(h(1),-h(2));
vr=dot(r,v)/norm(r); vt=norm(cross(r,v))/norm(r);
bins=[floor(node.t/43200),round(norm(r)/5000),round(vr/.5),round(vt/.5), ...
 round(inc/10),round(raan/30),floor(sum(node.visited)/7),kind,target,floor(log2(max(1,dt)))];
key=['policy|' sprintf('%d/',bins)];
end
