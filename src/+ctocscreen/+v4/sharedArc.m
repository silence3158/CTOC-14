function [children,report]=sharedArc(node,parent,eph,c,budget)
%SHAREDARC Joint two-encounter correction of the newly appended impulse.
% Same J2 STM/Newton equations as absorb; multi-start ranking and parent
% time boundary from the 2026-09-27 diagnostic. Partial search check only.
clock=tic; children={}; q=node.q; m=eph.model; M=numel(q.tau);
report=struct('status','no_candidate','events',0,'time_candidates',0,'near_events',0, ...
 'screened',0,'eligible',0,'attempts',0,'accepted',0,'far_accepted',0,'seconds',0, ...
 'failures',{{}},'trials',{{}});
if M~=numel(parent.q.tau)+1||budget<=0||q.tau(end)<parent.q.T, return; end
active=node.actual.distance_km<=1|parent.actual.distance_km<=1;
j=find(node.actual.distance_km<=1&abs(node.actual.witness_times_s-q.T)<1e-6,1);
if isempty(j), report.status='no_endpoint_anchor'; return; end
tau0=q.tau(end); u0=q.u(end,:).'; thetaJ=node.actual.witness_times_s(j);
if thetaJ<=tau0, report.status='no_positive_arc'; return; end
[ev,~]=ctocscreen.v4.events(node,eph,c); report.events=numel(ev);
latest=m.horizon_s-max(0,sum(~active)-2)*c.min_leg_s;
ev=ev([ev.time]<=latest); report.time_candidates=numel(ev);
report.near_events=sum([ev.miss_km]<=c.absorb_miss_km);
if isempty(ev)||toc(clock)>=budget, report.seconds=toc(clock); return; end
try
 xp=ctocscreen.v4.stateAt(node.trace,tau0,'post');
 [~,~,sol]=ctocscreen.v3Arc(xp,tau0,max([ev.time]),m,c,true);
catch err
 report.status='prescreen_propagation'; report.failures={err.identifier}; report.seconds=toc(clock); return
end
% Keep distinct absolute-time opportunities, including repeated targets.
rows=nan(numel(ev),5); scale=[.1;.1;.1;3600;3600;3600];
for k=1:numel(ev)
 if toc(clock)>=.65*budget, break; end
 [A,F]=equations(sol,thetaJ,ev(k).id,ev(k).time,u0);
 cond=rcond(A.*scale.'); report.screened=report.screened+1;
 if ~isfinite(cond)||cond<c.shared_rcond_min, continue; end
 dx=-((A.*scale.')\F).*scale;
 rows(k,:)=[norm(u0+dx(1:3))-norm(u0),max(abs(dx./scale)),norm(A*dx+F),cond,norm(dx(1:3))];
end
valid=find(all(isfinite(rows),2)&rows(:,1)<=c.absorb_max_extra_km_s&rows(:,3)<=c.absorb_residual_km);
report.eligible=numel(valid);
if isempty(valid), report.status='prescreen_rejected'; report.seconds=toc(clock); return; end
[~,s1]=sort(rows(valid,2)); [~,s2]=sort(rows(valid,1));
order=unique(reshape([valid(s1).';valid(s2).'],1,[]),'stable');
order=order(1:min(numel(order),c.shared_refine_candidates)); keys={};
for ik=order
 if toc(clock)>=budget-.02, break; end
 [qn,nr]=refine(ev(ik)); report.attempts=report.attempts+1;
 trial=struct('target',ev(ik).id,'event_s',ev(ik).time,'miss_km',ev(ik).miss_km, ...
  'predicted_extra',rows(ik,1),'scaled_step',rows(ik,2),'du_norm',rows(ik,5), ...
  'newton',nr,'accepted',false,'added',[],'actual_extra',NaN);
 if ~isempty(qn)
  [a,tr]=ctocscreen.v4.replay(qn,eph,c,false,node);
  n=numel(parent.q.tau); p=parent.q; qq=tr.q;
  prefix=numel(qq.tau)==n+1&&isequal(qq.x0,p.x0)&&isequal(qq.tau(1:n),p.tau) ...
   &&isequal(qq.u(1:n,:),p.u)&&qq.tau(end)>=p.T;
  added=find(a.distance_km<=c.search_radius_km&~active);
  if prefix&&a.initial_passed&&a.height_passed&&~strcmp(a.status,'propagation_failure') ...
    &&all(a.distance_km(active)<=1)&&a.distance_km(ev(ik).id)<=c.search_radius_km&&~isempty(added)
   key=ctocscreen.v4.controlKey(qq);
   if ~any(strcmp(keys,key))
    child=node; child.q=qq; child.q.witness=a.witness_times_s; child.actual=a; child.trace=tr;
    child.origin='shared_arc'; child.generation=parent.generation+1;
    child.attempts=0; child.zero_gain=0; child.heuristic_H=NaN;
    children{end+1}=child; keys{end+1}=key; %#ok<AGROW>
    trial.accepted=true; trial.added=added.'; trial.actual_extra=a.total_dv_km_s-node.actual.total_dv_km_s;
    report.accepted=report.accepted+1; report.far_accepted=report.far_accepted+(ev(ik).miss_km>c.absorb_miss_km);
   end
  else
   trial.replay_status=a.status; trial.prefix_valid=prefix;
  end
 end
 report.trials{end+1}=trial;
 if numel(children)>=c.shared_children, break; end
end
report.status='no_corrected_child'; if ~isempty(children), report.status='shared'; end
report.seconds=toc(clock);

 function [A,F]=equations(solution,tj,target,tn,u)
  [rj,vj]=ctocscreen.v3QueryTargets(eph,j,tj); [rn,vn]=ctocscreen.v3QueryTargets(eph,target,tn);
  yj=deval(solution,tj); yn=deval(solution,tn); Pj=reshape(yj(7:end),6,6); Pn=reshape(yn(7:end),6,6);
  F=[yj(1:3)-rj.';yn(1:3)-rn.'];
  A=[Pj(1:3,4:6),-Pj(1:3,1:3)*u,yj(4:6)-vj.',zeros(3,1); ...
     Pn(1:3,4:6),-Pn(1:3,1:3)*u,zeros(3,1),yn(4:6)-vn.'];
 end
 function [qn,nr]=refine(event)
  qn=[]; u=u0; tau=tau0; tj=thetaJ; tn=event.time; lo=parent.q.T;
  xpre=ctocscreen.v4.stateAt(node.trace,tau0,'pre');
  nr=struct('status','iteration_limit','iterations',0,'miss_km',[NaN NaN]);
  for it=1:c.absorb_newton_iterations
   nr.iterations=it;
   if toc(clock)>=budget, nr.status='budget'; return; end
   try
    if tau<=tau0, xm=ctocscreen.v4.stateAt(node.trace,tau,'pre'); else, xm=ctocscreen.v3Arc(xpre,tau0,tau,m,c); end
    [~,~,solution]=ctocscreen.v3Arc(xm+[0;0;0;u],tau,tn,m,c,true);
    [A,F]=equations(solution,tj,event.id,tn,u);
   catch err
    nr.status='propagation'; nr.failure=err.identifier; return
   end
   nr.miss_km=[norm(F(1:3)) norm(F(4:6))];
   if max(nr.miss_km)<=c.absorb_newton_km
    qn=q; qn.u(end,:)=u.'; qn.tau(end)=tau; qn.T=tn; qn.witness(j)=tj; qn.witness(event.id)=tn;
    nr.status='converged'; return
   end
   scaled=A.*scale.';
   if any(~isfinite(scaled(:)))||rcond(scaled)<c.shared_rcond_min, nr.status='ill_conditioned'; return; end
   dx=-(scaled\F).*scale;
   if any(~isfinite(dx)), nr.status='nonfinite_step'; return; end
   dx=dx*min(1,min(.2/max(norm(dx(1:3)),eps),1800/max(abs(dx(4:6)))));
   u=u+dx(1:3); tau=min(max(tau+dx(4),lo),tj-1);
   tj=min(m.horizon_s,max(tau+1,tj+dx(5))); tn=min(m.horizon_s,max(tj,tn+dx(6)));
  end
 end
end
