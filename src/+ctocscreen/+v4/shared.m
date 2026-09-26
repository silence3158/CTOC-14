function [child,report,warm]=shared(parent,eph,c,budget,scope)
%SHARED Jointly insert one or two encounters on an existing natural trajectory.
if nargin<5, scope='local'; end
child=[]; warm=[]; report=struct('status','no_opportunity','seconds',0,'iterations',0);
clock=tic; opp=ctocscreen.v4.opportunities(parent,eph,c);
if isempty(opp), report.seconds=toc(clock); return; end
active=find(parent.actual.distance_km<=1); theta=parent.actual.witness_times_s(active);
selected=opp(1); keep=1;
for k=2:numel(opp)
 if opp(k).id~=selected(1).id&&opp(k).arc==selected(1).arc
  selected(2)=opp(k); keep=2; break
 end
end
ids=[active;reshape([selected.id],[],1)]; theta=[theta;reshape([selected.time],[],1)];
seed=parent; seed.q.witness(ids)=theta;
if toc(clock)>=budget, report.status='budget_not_started'; report.seconds=toc(clock); return; end
[child,report,warm]=ctocscreen.v4.joint(seed,ids,theta,scope,selected(1).time,eph,c,budget-toc(clock));
report.added_target_ids=[selected.id]; report.shared_target_count=keep;
report.seconds=toc(clock);
if ~isempty(child), child.origin='shared_arc'; child.attempts=0; child.generation=parent.generation+1; end
end
