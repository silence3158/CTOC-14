function [child,report,warm,tabu]=shared(parent,eph,c,seconds,tabu)
%SHARED Jointly insert one or two encounters on an existing natural trajectory.
% A task is the control prefix up to the new encounter plus the added targets.
% Each executed failure escalates its scope (local, expanded, full); after
% shared_attempts failures, or one success, the task is skipped. This is a budget rule,
% not a proof that the encounter combination is infeasible.
if nargin<5||isempty(tabu), tabu=struct('keys',{{}},'failures',[]); end
child=[]; warm=[]; report=struct('status','no_opportunity','seconds',0,'iterations',0,'skipped_tasks',0);
clock=tic;
% Only the last arc: its encounters have no passive suffix, so opening the
% last burn gives as many free variables as new requirements (A3 design).
opp=ctocscreen.v4.opportunities(parent,eph,c,[],numel(parent.trace.arcs));
if isempty(opp), report.seconds=toc(clock); return; end
active=find(parent.actual.distance_km<=1); theta=parent.actual.witness_times_s(active);
levels={'local','expanded','full'}; chosen=[];
for first=1:numel(opp)
 selected=opp(first);
 for k=1:numel(opp)
  if k~=first&&opp(k).id~=selected(1).id&&opp(k).arc==selected(1).arc
   selected(2)=opp(k); break
  end
 end
 key=taskKey(parent.q,selected);
 index=find(strcmp(tabu.keys,key),1); failures=0;
 if ~isempty(index), failures=tabu.failures(index); end
 if failures<c.shared_attempts, chosen=selected; break; end
 report.skipped_tasks=report.skipped_tasks+1;
end
if isempty(chosen), report.status='all_tasks_exhausted'; report.seconds=toc(clock); return; end
selected=chosen; scope=levels{min(3,1+failures)}; allowance=seconds(min(3,1+failures));
ids=[active;reshape([selected.id],[],1)]; theta=[theta;reshape([selected.time],[],1)];
seed=parent; seed.q.witness(ids)=theta;
if toc(clock)>=allowance, report.status='budget_not_started'; report.seconds=toc(clock); return; end
[child,jr,warm]=ctocscreen.v4.joint(seed,ids,theta,scope,selected(1).time,eph,c,allowance-toc(clock));
jr.skipped_tasks=report.skipped_tasks; report=jr;
report.added_target_ids=[selected.id]; report.shared_target_count=numel(selected);
report.task_failures_before=failures; report.seconds=toc(clock);
% Only executed calls consume an attempt; resumable interruptions are queued
% elsewhere. A success on this prefix already produced its child, so the task
% is closed rather than repeated.
if report.iterations>0
 increment=1; if report.active_passed, increment=c.shared_attempts; end
 if isempty(index), tabu.keys{end+1}=key; tabu.failures(end+1)=min(c.shared_attempts,increment);
 else, tabu.failures(index)=min(c.shared_attempts,tabu.failures(index)+increment); end
end
if ~isempty(child), child.origin='shared_arc'; child.attempts=0; child.generation=parent.generation+1; child.heuristic_H=NaN; end
end
function key=taskKey(q,selected)
% Controls after the encounter do not define the insertion task.
keep=q.tau<=selected(1).time;
key=sprintf('%.17g,',[q.x0(:);q.tau(keep);reshape(q.u(keep,:),[],1);sort([selected.id]).']);
end
