function out=v3ScanArc(sol,targetQuery,c,witness,activeIds)
%V3SCANARC Search all targets; positive witnesses, no absence guarantee.
if nargin<4, witness=nan(35,1); end
if nargin<5, activeIds=1:35; end
out.distance_km=inf(35,1); out.time_s=nan(35,1); out.visited=false(35,1);
out.scope='requested targets searched; positive witnesses only, no proof of absence';
if isempty(activeIds), return; end
activeIds=activeIds(:).';
lo=sol.x(1); hi=sol.x(end);
times=linspace(lo,hi,max(3,ceil((hi-lo)/c.scan_step_s)+1));
x=deval(sol,times); [rt,~]=targetQuery(activeIds,times,'grid');
dist=reshape(vecnorm(rt-permute(x(1:3,:),[3 1 2]),2,2),numel(activeIds),[]);
for index=1:numel(activeIds)
 id=activeIds(index); row=dist(index,:); [d,ix]=min(row); tt=times(ix);
 local=find(row(2:end-1)<=row(1:end-2)&row(2:end-1)<=row(3:end))+1;
 local=unique([ix local]);
 for j=local
  a=times(max(1,j-1)); b=times(min(numel(times),j+1));
  [tq,dq]=fminbnd(@distance,a,b,optimset('TolX',1e-6,'Display','off'));
  if dq<d, d=dq; tt=tq; end
 end
 if isfinite(witness(id))&&witness(id)>=lo&&witness(id)<=hi
  dw=distance(witness(id)); if dw<d, d=dw; tt=witness(id); end
 end
 out.distance_km(id)=d; out.time_s(id)=tt;
end
out.visited=out.distance_km<=1;
 function d=distance(t)
  y=deval(sol,t); r=targetQuery(id,t,'pairs'); d=norm(y(1:3).'-r);
 end
end
