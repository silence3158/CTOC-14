function [child,report]=absorb(node,eph,c,budget,tabu)
%ABSORB Fold a following encounter into the node's last burns (one burn, two targets).
% After a transfer arrives at target j, remaining targets that cross the
% inspector plane soon afterwards with a small miss are candidates. A linear
% prescreen (MATH_SPEC eq. 29; eq. 19 for the burn-time column) solves
%   e_j + S_j*d + w_j*dtheta_j = 0,  e_n + S_n*d + w_n*dtheta_n = 0,
% d = [du; dtau] for the last burn, and predicts the new pulse norm. Only a
% predicted cost increase below absorb_max_extra_km_s is kept; it is
% then refined by Newton-Raphson on the real propagation (refine below).
% The prescreen is a ranking heuristic; fixed-control replay decides acceptance.
if nargin<5, tabu={}; end
clock=tic; child=[]; m=eph.model; q=node.q;
report=struct('status','no_candidate','candidates',0,'screened',0,'target',NaN,'predicted_extra',NaN, ...
 'prescreen_residual_km',NaN,'seconds',0,'joint',[],'added',[],'newton',[]);
M=numel(q.tau); if M==0||toc(clock)>=budget, return; end
[ev,info]=ctocscreen.v4.events(node,eph,c);
if isempty(ev), report.seconds=toc(clock); return; end
left=sum(node.actual.distance_km>1); latest=m.horizon_s-max(0,left-2)*c.min_leg_s;
ok=[ev.miss_km]<=c.absorb_miss_km&[ev.time]<=latest; ev=ev(ok);
[~,o]=sort([ev.miss_km]); ev=ev(o); seen=false(35,1); keep=false(size(ev));
for k=1:numel(ev)
 key=sprintf('%s|%d',ctocscreen.v4.controlKey(q),ev(k).id);
 if ~seen(ev(k).id)&&~any(strcmp(tabu,key)), keep(k)=true; seen(ev(k).id)=true; end
end
ev=ev(keep); ev=ev(1:min(end,c.absorb_candidates)); report.candidates=numel(ev);
if isempty(ev), report.seconds=toc(clock); return; end
% Linear prescreen along the last arc with one variational integration.
tauM=q.tau(M); uM=q.u(M,:).'; xplus=ctocscreen.v4.stateAt(node.trace,tauM,'post');
j=find(node.actual.distance_km<=1&abs(node.actual.witness_times_s-q.T)<1e-6,1);
thetaJ=q.T; if isempty(j), [~,j]=max(node.actual.witness_times_s.*(node.actual.distance_km<=1)); thetaJ=node.actual.witness_times_s(j); end
tEnd=max([ev.time]);
try
 [~,~,sol]=ctocscreen.v3Arc(xplus,tauM,tEnd,m,c,true);
catch err
 report.status='prescreen_failure'; report.failure_id=err.identifier; report.seconds=toc(clock); return
end
best=[]; bestExtra=Inf;
for k=1:numel(ev)
 [row,res]=screen(ev(k)); report.screened=report.screened+1;
 if res<=c.absorb_residual_km&&row.extra<bestExtra, best=row; bestExtra=row.extra; end
end
if isempty(best)||bestExtra>c.absorb_max_extra_km_s
 report.status='prescreen_rejected'; report.predicted_extra=bestExtra; report.seconds=toc(clock); return
end
report.target=best.id; report.predicted_extra=best.extra; report.prescreen_residual_km=best.res;
% Newton-Raphson on the square two-encounter system (Xia et al. 2022 solve
% their one-impulse two-target reduction the same way): unknowns du, tau,
% thetaJ, thetaN; residual both encounter miss vectors; real J2 propagation.
[qn,newton]=refine(best);
report.newton=newton;
if isempty(qn), report.status=['newton_',newton.status]; report.seconds=toc(clock); return; end
[a,tr]=ctocscreen.v4.replay(qn,eph,c,false,node);
active=find(node.actual.distance_km<=1);
if strcmp(a.status,'propagation_failure')||~a.height_passed||~a.initial_passed ...
  ||any(a.distance_km(active)>c.search_radius_km)||a.distance_km(best.id)>c.search_radius_km
 report.status='replay_rejected'; report.seconds=toc(clock); return
end
child=node; child.q=tr.q; child.q.witness=a.witness_times_s; child.actual=a; child.trace=tr;
child.origin='absorbed_encounter'; child.attempts=0; child.zero_gain=0; child.heuristic_H=NaN;
child.generation=node.generation+1; report.added=find(a.distance_km<=1&node.actual.distance_km>1).';
report.status='absorbed'; report.seconds=toc(clock);
 function [qn,info]=refine(row)
  info=struct('status','iteration_limit','iterations',0,'miss_km',[NaN NaN]);
  qn=[]; u=uM; tau=tauM; tJ=thetaJ; tN=row.time;
  lo=0; if M>1, lo=q.tau(M-1); end
  xpre=ctocscreen.v4.stateAt(node.trace,tauM,'pre');
  for it=1:c.absorb_newton_iterations
   info.iterations=it;
   if toc(clock)>=budget, info.status='budget'; return; end
   try
    if tau<=tauM, xm=ctocscreen.v4.stateAt(node.trace,tau,'pre'); else, xm=ctocscreen.v3Arc(xpre,tauM,tau,m,c); end
    [~,~,s2]=ctocscreen.v3Arc(xm+[0;0;0;u],tau,max(tJ,tN),m,c,true);
   catch
    info.status='propagation'; return
   end
   [rj,vj]=ctocscreen.v3QueryTargets(eph,j,tJ); [rn,vn]=ctocscreen.v3QueryTargets(eph,row.id,tN);
   yj=deval(s2,tJ); yn=deval(s2,tN); Pj=reshape(yj(7:end),6,6); Pn=reshape(yn(7:end),6,6);
   F=[yj(1:3)-rj.';yn(1:3)-rn.']; info.miss_km=[norm(F(1:3)) norm(F(4:6))];
   if max(info.miss_km)<=c.absorb_newton_km
    qn=q; qn.u(M,:)=u.'; qn.tau(M)=tau; qn.T=tN; qn.witness(j)=tJ; qn.witness(row.id)=tN;
    info.status='converged'; return
   end
   A=[Pj(1:3,4:6),-Pj(1:3,1:3)*u,yj(4:6)-vj.',zeros(3,1);Pn(1:3,4:6),-Pn(1:3,1:3)*u,zeros(3,1),yn(4:6)-vn.'];
   sc=[.1*ones(3,1);3600;3600;3600]; dx=-((A.*sc.')\F).*sc;
   % Damp to a bounded step (0.2 km/s, 30 min) and keep the event order.
   dx=dx*min(1,min(.2/max(norm(dx(1:3)),eps),1800/max(abs(dx(4:6)))));
   u=u+dx(1:3); tau=min(max(tau+dx(4),lo),tJ-1); tJ=max(tau+1,tJ+dx(5)); tN=min(m.horizon_s,max(tJ,tN+dx(6)));
  end
 end
 function [row,res]=screen(e)
  % Rows: 3 for target j, 3 for target n; columns du(3), dtau, dthetaJ, dthetaN.
  [rj,vj]=ctocscreen.v3QueryTargets(eph,j,thetaJ); [rn,vn]=ctocscreen.v3QueryTargets(eph,e.id,e.time);
  yj=deval(sol,thetaJ); yn=deval(sol,e.time); Pj=reshape(yj(7:end),6,6); Pn=reshape(yn(7:end),6,6);
  aj=ctocscreen.v3Force(thetaJ,yj(1:3),m); an=ctocscreen.v3Force(e.time,yn(1:3),m); %#ok<NASGU>
  % d r(theta)/d tau = -Phi_rr * u (eq. 19 with f^- - f^+ = [-u; 0]).
  A=[Pj(1:3,4:6),-Pj(1:3,1:3)*uM,yj(4:6)-vj.',zeros(3,1); ...
     Pn(1:3,4:6),-Pn(1:3,1:3)*uM,zeros(3,1),yn(4:6)-vn.'];
  b=-[yj(1:3)-rj.';yn(1:3)-rn.'];
  s=[ones(3,1)*.1;3600;3600;3600]; x=(A.*s.')\b; x=x.*s; res=norm(A*x-b);
  row=struct('id',e.id,'time',e.time,'extra',norm(uM+x(1:3))-norm(uM),'res',res,'dtau',x(4),'du',norm(x(1:3)));
 end
end
