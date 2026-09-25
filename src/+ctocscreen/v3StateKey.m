function key=v3StateKey(n,kind,target,dt,c)
%V3STATEKEY Configurable bins of real Cartesian state and exact visited set.
if nargin<5, c=ctocscreen.v3Defaults(); end
s=c.pheromone_state_scales; x=n.state(:).';
bins=[floor(n.t/s(1)),round(x(1:3)/s(2)),round(x(4:6)/s(3)), ...
 kind,target,floor(log2(max(1,dt)))];
key=[sprintf('%d/',bins),char('0'+n.visited(:).')];
end
