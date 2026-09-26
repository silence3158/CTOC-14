function opp=opportunities(node,eph,c,ids,arcs)
%OPPORTUNITIES Retain multiple local windows, including different revolutions.
% arcs restricts the scan (default: all arcs of the trace).
if nargin<4||isempty(ids), ids=find(node.actual.distance_km>1); end
if nargin<5||isempty(arcs), arcs=1:numel(node.trace.arcs); end
opp=struct('id',{},'time',{},'distance',{},'arc',{});
if isempty(ids), return; end
for a=arcs(:).'
 sol=node.trace.arcs{a}; lo=sol.x(1); hi=sol.x(end);
 tt=linspace(lo,hi,max(3,ceil((hi-lo)/c.scan_step_s)+1)); xx=deval(sol,tt);
 rt=ctocscreen.v3QueryTargets(eph,ids,tt,'grid');
 dd=reshape(vecnorm(rt-permute(xx(1:3,:),[3 1 2]),2,2),numel(ids),[]);
 for j=1:numel(ids)
  [~,best]=min(dd(j,:)); minima=find(dd(j,2:end-1)<=dd(j,1:end-2)&dd(j,2:end-1)<=dd(j,3:end))+1;
  for ix=unique([best,minima])
   if dd(j,ix)>2*c.discovery_radius_km, continue; end
   ta=tt(max(1,ix-1)); tb=tt(min(numel(tt),ix+1));
   [t,d]=fminbnd(@miss,ta,tb,optimset('TolX',1e-4,'Display','off'));
   if d<=c.discovery_radius_km
    opp(end+1)=struct('id',ids(j),'time',t,'distance',d,'arc',a); %#ok<AGROW>
   end
  end
 end
end
if ~isempty(opp), [~,order]=sort([opp.distance]); opp=opp(order); end
 function d=miss(t)
  yy=deval(sol,t); rr=ctocscreen.v3QueryTargets(eph,ids(j),t);
  d=norm(yy(1:3)-rr.');
 end
end
