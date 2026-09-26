function [child,report,warm]=restructure(parent,eph,c,index,budget)
%RESTRUCTURE Changed topology is a new task; actual old coverage is re-evaluated.
% index counts executed structure services, so the five operators rotate
% independently of root or full-history rounds.
clock=tic; child=[]; warm=[]; report=struct('status','no_structure','seconds',0,'iterations',0);
q=parent.q; M=numel(q.tau); if q.T<=0, return; end
kind=mod(index-1,5);
if M==0, kind=1; end
if kind==4
 [child,report]=ctocscreen.v4.rebuild(parent,eph,c,ceil(index/5));
 return
end
ids=find(parent.actual.distance_km<=1); theta=parent.actual.witness_times_s(ids);
switch kind
 case 0
  [~,k]=max(vecnorm(q.u,2,2)); focus=q.tau(k); q.tau(k)=[]; q.u(k,:)=[]; name='delete_expensive';
 case 1
  if M>=c.max_maneuvers, return; end
  edges=[0;q.tau;q.T]; [~,k]=max(diff(edges)); focus=(edges(k)+edges(k+1))/2;
  q.tau(end+1,1)=focus; q.u(end+1,:)=zeros(1,3); name='insert_zero';
 case 2
  k=1+mod(index-1,M); focus=q.tau(k);
  q.tau(k)=max(0,min(q.T-1,focus+.2*c.window_radius_s)); name='move_event';
 otherwise
  [~,k]=max(vecnorm(q.u,2,2)); focus=q.tau(k);
  if k<M
   q.u(k,:)=q.u(k,:)+q.u(k+1,:); q.u(k+1,:)=[]; q.tau(k+1)=[]; name='merge_arcs';
  else
   q.T=min(eph.model.horizon_s,q.T+c.lookahead_s/2); name='alternate_windows';
  end
end
[q.tau,order]=sort(q.tau); q.u=q.u(order,:);
[a,tr]=ctocscreen.v4.replay(q,eph,c);
if strcmp(a.status,'propagation_failure'), report.status=a.status; report.seconds=toc(clock); return; end
seed=parent; seed.q=q; seed.actual=a; seed.trace=tr;
if strcmp(name,'alternate_windows')
 opp=ctocscreen.v4.opportunities(seed,eph,c,ids);
 if ~isempty(opp)
  j=find(ids==opp(1).id,1); theta(j)=opp(1).time;
 end
end
if toc(clock)>=budget, report.status='budget_not_started'; report.seconds=toc(clock); return; end
[child,report,warm]=ctocscreen.v4.joint(seed,ids,theta,'full',focus,eph,c,budget-toc(clock));
report.action=name; report.seconds=toc(clock);
if ~isempty(child), child.origin=name; child.attempts=0; child.generation=parent.generation+1; child.heuristic_H=NaN; end
end
