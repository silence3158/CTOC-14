function [ev,info]=events(node,eph,c)
%EVENTS Encounter proposals where remaining targets meet the inspector plane.
% Ballistic coast from the node end over a finite window (search approximation).
% Crossing events: sign changes of hhat(t)'*r_j(t). Near-coplanar targets use
% in-plane radius matches instead. These are proposals, never visit evidence.
% Idea after the traversing-point method (Zhang & Zhou 2013, local PDF);
% timing and J2 are handled by real propagation, not by a fitted ellipse.
m=eph.model; q=node.q; T=q.T; mu=m.mu;
ev=struct('id',{},'time',{},'radius',{},'miss_km',{},'dv_est',{},'kind',{});
% H defaults to a large value: a node without events is not cheap to finish.
info=struct('sol',[],'t_start',T,'t_end',T,'H',NaN,'target_estimate',nan(35,1),'seconds',0,'failure_id','');
clock=tic; remaining=find(node.actual.distance_km>1);
W=min(c.event_window_s,m.horizon_s-T);
if isempty(remaining), info.H=0; info.seconds=toc(clock); return; end
info.H=c.no_event_cost_km_s*numel(remaining);
if W<=c.min_flight_s, info.seconds=toc(clock); return; end
% Integrate from the last real impulse, as fixed-control replay does. The
% look-ahead may hit the Earth; that only means no ballistic events here.
if isempty(q.tau), a0=0; xa=q.x0; else, a0=q.tau(end); xa=ctocscreen.v4.stateAt(node.trace,a0,'post'); end
try
 [~,~,sol]=ctocscreen.v3Arc(xa,a0,T+W,m,c,false,true);
catch err
 info.failure_id=err.identifier; info.seconds=toc(clock); return
end
% A look-ahead that dips below 200 km is truncated at the first crossing.
tg=unique([T:c.event_grid_s:T+W,T+W]); X=deval(sol,tg); low=find(vecnorm(X(1:3,:))<m.re+200,1);
if ~isempty(low)
 if low<=2, info.seconds=toc(clock); return; end
 tg=tg(1:low-1); W=tg(end)-T;
end
nt=numel(tg); ni=numel(remaining);
X=X(:,1:nt); r=X(1:3,:); v=X(4:6,:); hr=cross(r,v); hn=hr./vecnorm(hr);
[rt,vt]=ctocscreen.v3QueryTargets(eph,remaining,tg,'grid');
s=reshape(sum(rt.*permute(hn,[3 1 2]),2),ni,nt);
% Osculating conic radius toward each target direction, for coplanar matches.
R=vecnorm(r); V2=sum(v.^2,1); p=sum(hr.^2,1)/mu; e=((V2-mu./R).*r-sum(r.*v,1).*v)/mu;
ht=cross(rt(:,:,1),vt(:,:,1),2); ht=ht./vecnorm(ht,2,2);
planeAngle=acosd(max(-1,min(1,ht*hn(:,1))));
times=[]; ids=[]; kinds=[];
for j=1:ni
 if planeAngle(j)<c.coplanar_deg
  u=reshape(rt(j,:,:),3,nt); Rj=vecnorm(u); u=u./Rj; u=u-sum(u.*hn,1).*hn; u=u./vecnorm(u);
  f=p./(1+sum(e.*u,1))-Rj; kind=2;
 else
  f=s(j,:); kind=1;
 end
 k=find(f(1:end-1).*f(2:end)<=0&(f(1:end-1)~=0|f(2:end)~=0));
 w=f(k)./(f(k)-f(k+1)); t=tg(k)+w.*(tg(k+1)-tg(k));
 times=[times,t]; ids=[ids,repmat(remaining(j),1,numel(t))]; kinds=[kinds,repmat(kind,1,numel(t))]; %#ok<AGROW>
end
keep=times>=T+c.min_flight_s; times=times(keep); ids=ids(keep); kinds=kinds(keep);
info.sol=sol; info.t_end=T+W;
if isempty(times), info.seconds=toc(clock); return; end
Xe=deval(sol,times); [re,~]=ctocscreen.v3QueryTargets(eph,ids,times);
miss=vecnorm(Xe(1:3,:).'-re,2,2);
% One tangential burn opposite the target direction, apse-to-apse (Hohmann-type).
% A heuristic for ranking only: it ignores phasing and is no lower bound.
x=X(:,1); r0=x(1:3); v0=x(4:6); R0=norm(r0); a=1/(2/R0-dot(v0,v0)/mu);
h0=cross(r0,v0); p0=dot(h0,h0)/mu; h0=h0/norm(h0); e0=((dot(v0,v0)-mu/R0)*r0-dot(r0,v0)*v0)/mu;
Rj=vecnorm(re,2,2); u=re./Rj; off=u*h0; up=u-off*h0.'; up=up./vecnorm(up,2,2);
ro=p0./(1-up*e0); vo=sqrt(max(0,mu*(2./ro-1/a))); vn=sqrt(max(0,mu*(2./ro-2./(ro+Rj))));
dv=abs(vn-vo)+vo.*abs(asin(max(-1,min(1,off))));
for k=1:numel(times)
 ev(end+1)=struct('id',ids(k),'time',times(k),'radius',Rj(k),'miss_km',miss(k), ...
  'dv_est',dv(k),'kind',kinds(k)); %#ok<AGROW>
end
best=inf(35,1);
for k=1:numel(ev), best(ev(k).id)=min(best(ev(k).id),ev(k).dv_est); end
est=best(remaining); finite=est(isfinite(est));
if ~isempty(finite), est(~isfinite(est))=max(finite); end
info.target_estimate(remaining)=est;
if ~isempty(finite)
 sorted=sort(est); n=max(1,ceil(c.heuristic_fraction*numel(sorted)));
 info.H=c.heuristic_weight*numel(remaining)*mean(sorted(1:n));
end
info.seconds=toc(clock);
end
