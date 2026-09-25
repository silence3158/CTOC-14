function [pairs,report]=v3TimeCandidates(node,eph,c,stream,pheromone,clock,coverage)
%V3TIMECANDIDATES Rank target/time pairs by J2 position-only STM corrections.
if nargin<7, coverage=false; end
if c.cost_guidance_enabled
 [pairs,report]=ctocscreen.v3CostCandidates(node,eph,c,stream,pheromone,clock); return;
end
ids=find(~node.visited); pairs={};
report=struct('propagations',0,'pairs_evaluated',0,'failures',{{}},'method','j2_stm_position_only', ...
 'pheromone_queries',0,'nonneutral_queries',0,'prefix_guidance_hits',0);
if isempty(ids), return; end
remaining=eph.model.horizon_s-node.t;
grid=c.action_durations_s(:)*(0.85+0.3*rand(stream));
if coverage
 % Reserve time for the remaining visits as a construction heuristic only.
 slot=remaining/numel(ids);
 grid=min(slot,c.cold_durations_s(:));
end
grid=unique(min(remaining,max(c.action_min_duration_s,grid)));
if numel(grid)>c.time_grid_count
 grid=grid(unique(round(linspace(1,numel(grid),c.time_grid_count))));
end
if node.t==0&&isfield(node,'seed_target')&&node.seed_target>0&&isfinite(node.seed_duration_s)
 grid=unique([grid;min(remaining,node.seed_duration_s)]);
end
if isfield(node,'guidance_times_s')
 hints=node.guidance_times_s(ids)-node.t; hints=hints(isfinite(hints)&hints>0&hints<=remaining);
 if numel(hints)>c.time_grid_count, hints=hints(unique(round(linspace(1,numel(hints),c.time_grid_count)))); end
 grid=unique([grid;hints(:)]);
end
grid=grid(grid>=c.action_min_duration_s & grid<=remaining);
for dt=grid(:).'
 if toc(clock)>=c.budget_s, break; end
 try
  [y,P]=ctocscreen.v3Arc(node.state,node.t,node.t+dt,eph.model,c,true);
  goals=ctocscreen.v3QueryTargets(eph,ids,node.t+dt);
  B=P(1:3,4:6); miss=goals.'-y(1:3);
  corrections=[B;1e-5*eye(3)]\[miss;zeros(3,numel(ids))];
  score=vecnorm(corrections,2,1); report.propagations=report.propagations+1;
  for k=1:numel(ids)
   if ~isfinite(score(k)), continue; end
   key=ctocscreen.v3StateKey(node,3,ids(k),dt,c);
   [ph,~,feedback]=ctocscreen.v3Pheromone('get',pheromone,key,0,c);
   report.pheromone_queries=report.pheromone_queries+1;
   report.nonneutral_queries=report.nonneutral_queries+(ph>c.pheromone_baseline);
   report.prefix_guidance_hits=report.prefix_guidance_hits+(feedback.prefix_increment>0);
   pairs{end+1}=struct('target',ids(k),'dt',dt,'estimate',score(k),'pheromone',ph, ...
    'natural_miss_km',norm(miss(:,k)),'weight',0);
  end
 catch err
  report.failures{end+1}=struct('id',err.identifier,'message',err.message,'dt',dt);
 end
end
report.pairs_evaluated=numel(pairs);
if isempty(pairs), return; end
scores=cellfun(@(p)p.estimate,pairs); [~,order]=sort(scores); ranks=zeros(size(order)); ranks(order)=0:numel(order)-1;
for k=1:numel(pairs)
 eta=(1-ranks(k)/numel(pairs))^c.time_rank_power;
 boost=1;
 if isfield(node,'priority_targets')&&ismember(pairs{k}.target,node.priority_targets), boost=2; end
 urgency=1;
 if coverage, urgency=min(1,slot/pairs{k}.dt)^2; end
 pairs{k}.weight=boost*urgency*pairs{k}.pheromone^c.pheromone_alpha*eta^c.heuristic_beta;
end
end
