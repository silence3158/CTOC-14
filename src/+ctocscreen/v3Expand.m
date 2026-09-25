function [children,failures,report]=v3Expand(node,eph,c,stream,pheromone,deadline)
%V3EXPAND Preserve multiple target/time/arrival-state successors.
children={}; failures={}; started=toc(deadline);
% Rotate the mandatory families using the recorded RNG, so short budgets do
% not always privilege coasts and random burns over guided construction.
order=randperm(stream,3);
report=struct('time_search',{{}},'generated',0,'retained',0, ...
 'action_order',order,'attempted_kinds',[],'guided_actions',0,'guided_processed',0, ...
 'reused_arcs',0,'elapsed_s',0,'budget_overrun_s',0,'policy_feedback',{{}});
for trial=1:c.actions_per_node
 if toc(deadline)>=c.budget_s, break; end
 try
  remaining=eph.model.horizon_s-node.t;
  dt=min(remaining,c.action_durations_s(randi(stream,numel(c.action_durations_s)))*(.7+.6*rand(stream)));
  if dt<c.action_min_duration_s, continue; end
  kind=order(min(trial,3)); wait=0;
  if trial>3
   wait=.8*dt*rand(stream); free=node;
   if wait>0&&pheromone.Count>0
    [x,~,sol]=ctocscreen.v3Arc(node.state,node.t,node.t+wait,eph.model,c);
    scan=ctocscreen.v3ScanArc(sol,@(ids,t,mode)ctocscreen.v3QueryTargets(eph,ids,t,mode),c,nan(35,1));
    free.state=x.'; free.t=node.t+wait; free.visited=node.visited|scan.visited;
   end
   weights=ones(1,3)*c.pheromone_baseline;
   weights(1)=actionWeight(node,1,0,dt);
   weights(2)=actionWeight(free,2,0,dt-wait);
   ids=find(~node.visited);
   if ~isempty(ids)
    weights(3)=mean(arrayfun(@(id)actionWeight(node,3,id,dt),ids));
   end
   weights=weights.^c.pheromone_alpha;
   weights=(1-c.exploration)*weights/sum(weights)+c.exploration/3;
   kind=find(cumsum(weights)>=rand(stream),1);
  end
  if toc(deadline)>=c.budget_s, break; end
  report.attempted_kinds(end+1)=kind;
  if kind==3
   local=c;
   local.budget_s=toc(deadline)+max(0,c.budget_s-toc(deadline))/(c.actions_per_node-trial+1);
   [actions,search]=ctocscreen.v3TimeSearch(node,eph,local,stream,pheromone,deadline);
   report.time_search{end+1}=search;
   report.guided_actions=report.guided_actions+numel(actions);
   [guided,failed,finished]=ctocscreen.v3GuidedChildren(node,actions,eph,c);
   report.guided_processed=report.guided_processed+finished.processed;
   report.reused_arcs=report.reused_arcs+finished.reused_arcs;
   report.generated=report.generated+finished.generated;
   failures=[failures failed]; children=[children guided];
  else
   dv=zeros(3,1);
   if kind==2
    if trial<=3, wait=.8*dt*rand(stream); end
    d=randn(stream,3,1); dv=d/norm(d)*c.random_pulse_scales(randi(stream,numel(c.random_pulse_scales)));
   end
   children{end+1}=ctocscreen.v3ApplyAction(node,node.t+wait,dv,node.t+dt,0,eph,c);
   a=children{end}.last_action;
   report.policy_feedback{end+1}=ctocscreen.v3PolicyFeedback(pheromone,a.origin,0,a.dt,a.dv,c);
   report.generated=report.generated+1;
  end
 catch err
  failures{end+1}=struct('id',err.identifier,'message',err.message);
 end
end
report.retained=numel(children);
report.elapsed_s=toc(deadline)-started;
report.budget_overrun_s=max(0,toc(deadline)-c.budget_s);
 function value=actionWeight(origin,kind,target,dt)
  exact=ctocscreen.v3Pheromone('get',pheromone,ctocscreen.v3StateKey(origin,kind,target,dt,c),0,c);
  coarse=ctocscreen.v3Pheromone('get',pheromone,ctocscreen.v3PolicyKey(origin,kind,target,dt),0,c);
  value=exact*coarse/c.pheromone_baseline;
 end
end
