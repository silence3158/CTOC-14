function [node,stats,failures,candidates]=v3Complete(input,eph,c,stream,pheromone,screenPrefixes,coverage)
%V3COMPLETE Bounded beam completion with multiple independently screened outputs.
if nargin<6, screenPrefixes=true; end
if nargin<7, coverage=false; end
if coverage, c.completion_beam_width=min(c.completion_beam_width,c.cold_beam_width); end
clock=tic; if iscell(input), roots=input; else, roots={input}; end
frontier=ctocscreen.v3SelectBeam(roots,c.completion_beam_width,c); original=frontier{1};
finished={}; stopped={}; failures={};
% Scope the cache to this immutable ephemeris/configuration and completion call.
checks=containers.Map('KeyType','char','ValueType','any');
stats=struct('steps',0,'elapsed_s',0,'start_visits',sum(original.visited),'end_visits',0, ...
 'input_count',numel(frontier),'layer_widths',[],'expansion_history',{{}},'candidate_checks',{{}}, ...
 'consistency_checks',0,'consistency_cache_hits',0,'prefix_check_seconds',0, ...
 'output_count',0,'prefix_screening',screenPrefixes,'coverage_completion',coverage);
for step=1:c.completion_steps
 if toc(clock)>=.8*c.completion_seconds||isempty(frontier), break; end
 next={};
 for k=1:numel(frontier)
  parent=frontier{k};
  if all(parent.visited), finished{end+1}=parent; continue; end
  if parent.t>=eph.model.horizon_s, stopped{end+1}=parent; continue; end
  if toc(clock)>=.8*c.completion_seconds
   stopped=[stopped frontier(k:end)]; break;
  end
  local=c; local.actions_per_node=max(3,c.actions_per_node);
  local.budget_s=toc(clock)+max(0,.8*c.completion_seconds-toc(clock))/(numel(frontier)-k+1);
  if coverage
   [children,failed,expansion]=ctocscreen.v3TailExpand(parent,eph,local,stream,pheromone,clock);
  else
   [children,failed,expansion]=ctocscreen.v3Expand(parent,eph,local,stream,pheromone,clock);
  end
  failures=[failures failed]; stats.expansion_history{end+1}=expansion;
  if isempty(children), stopped{end+1}=parent; else, next=[next children]; end
 end
 if isempty(next), frontier={}; break; end
 frontier=ctocscreen.v3SelectBeam(next,c.completion_beam_width,c);
 stats.steps=stats.steps+1; stats.layer_widths(end+1)=numel(frontier);
 if mod(step,c.consistency_stride)==0
  frontier=screenNodes(frontier);
 end
end
pool=ctocscreen.v3SelectBeam([finished frontier stopped],c.completion_outputs,c);
if isempty(pool), pool={original}; end
candidates=screenNodes(pool);
if isempty(candidates)
 node=original; stats.replay_status='rejected_new_prefix';
 stats.consistency=struct('passed',false,'truncated',false); stats.incremental_visits=sum(node.visited);
else
 coverage=cellfun(@(n)sum(n.visited),candidates); costs=cellfun(@(n)n.J,candidates);
 [~,order]=sortrows([-coverage(:),costs(:)]); candidates=candidates(order); node=candidates{1};
 stats.consistency=node.completion_check; stats.replay_status='passed_height';
 if ~screenPrefixes, stats.replay_status='proposal_pending_full_recovery'; end
 stats.incremental_visits=sum(node.visited);
end
stats.end_visits=sum(node.visited); stats.output_count=numel(candidates); stats.elapsed_s=toc(clock);
stats.budget_overrun_s=max(0,stats.elapsed_s-c.completion_seconds);

 function output=screenNodes(nodes)
  output={};
  for j=1:numel(nodes)
   n=nodes{j};
   if n.t<=0, continue; end
   if ~screenPrefixes
    n.completion_check=struct('passed',false,'truncated',false);
    output{end+1}=n; continue;
   end
   key=[ctocscreen.v3ScheduleKey(n.schedule) '|w=' mat2str(n.schedule.witness_times_s,17)];
   if isKey(checks,key)
    saved=checks(key); stable=saved.schedule; check=saved.screen; consistency=saved.report;
    stats.consistency_cache_hits=stats.consistency_cache_hits+1;
   else
    [stable,check,consistency]=ctocscreen.v3StablePrefix(n.schedule,eph,c);
    checks(key)=struct('schedule',stable,'screen',check,'report',consistency);
    stats.consistency_checks=stats.consistency_checks+1;
    stats.prefix_check_seconds=stats.prefix_check_seconds+consistency.elapsed_s;
   end
   stats.candidate_checks{end+1}=consistency;
   if consistency.passed
    n.schedule=stable; n.t=stable.duration_s; n.state=check.final_state;
    n.visited=check.distance_km<=1; n.J=check.total_dv_km_s; n.estimate=n.J;
    n.inclination_penalty=consistency.inclination_penalty;
    n.estimate=n.estimate+c.plane_penalty_km_s*n.inclination_penalty;
    if c.cost_guidance_enabled
     [radial,shape]=ctocscreen.v3ContinuationEstimate(n.state,n.visited,eph);
     n.estimate=n.estimate+c.radial_tour_weight*radial+c.arrival_shape_weight*shape;
    end
    n.completion_check=consistency; output{end+1}=n;
   else
    failures{end+1}=struct('id','ctocscreen:v3:completionReplay','message',check.status);
   end
  end
 end
end
