function [children,failures,report]=v3TailExpand(node,eph,c,stream,pheromone,clock)
%V3TAILEXPAND CEA-inspired tail addition; outer beam still explores free arcs.
started=toc(clock); local=c; local.guided_attempts=min(2,c.guided_attempts);
if c.cost_guidance_enabled, local.guided_attempts=max(4,c.guided_attempts); end
local.guided_branches=min(2,c.guided_branches); local.time_refine_iterations=0;
[actions,search]=ctocscreen.v3TimeSearch(node,eph,local,stream,pheromone,clock,true);
searches={search};
[children,failures,done]=ctocscreen.v3GuidedChildren(node,actions,eph,c);
children=children(cellfun(@(n)sum(n.visited)>sum(node.visited),children));
if isempty(children)&&toc(clock)<c.budget_s
 % Short transfers condition the cold construction; longer arcs remain available.
 [more,fallback]=ctocscreen.v3TimeSearch(node,eph,local,stream,pheromone,clock,false);
 searches{end+1}=fallback;
 [children,failed,extra]=ctocscreen.v3GuidedChildren(node,more,eph,c);
 children=children(cellfun(@(n)sum(n.visited)>sum(node.visited),children));
 actions=[actions more]; failures=[failures failed];
 done.generated=done.generated+extra.generated; done.processed=done.processed+extra.processed;
 done.reused_arcs=done.reused_arcs+extra.reused_arcs;
end
feedback={};
if c.cost_guidance_enabled&&toc(clock)<c.budget_s
 [transitions,failed]=ctocscreen.v3OrbitShapeChildren(node,eph,c);
 for j=1:numel(transitions)
  a=transitions{j}.last_action;
  feedback{end+1}=ctocscreen.v3PolicyFeedback(pheromone,a.origin,0,a.dt,a.dv,c);
 end
 children=[children transitions]; failures=[failures failed];
end
report=struct('time_search',{searches},'generated',done.generated,'retained',numel(children), ...
 'guided_actions',numel(actions),'guided_processed',done.processed,'reused_arcs',done.reused_arcs, ...
 'elapsed_s',toc(clock)-started,'budget_overrun_s',max(0,toc(clock)-c.budget_s), ...
 'method','coverage_tail_addition');
report.policy_feedback=feedback;
end
