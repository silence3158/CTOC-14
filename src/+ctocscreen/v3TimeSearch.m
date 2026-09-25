function [actions,report]=v3TimeSearch(node,eph,c,stream,pheromone,clock,coverage)
%V3TIMESEARCH Multi-time branching plus bounded pattern refinement of true DV.
if nargin<7, coverage=false; end
actions={}; report=struct('screen',[],'trials',{{}},'selected_pairs',{{}},'local_improvements',0,'policy_feedback',{{}});
screen=c; screen.budget_s=min(c.budget_s,toc(clock)+.25*max(0,c.budget_s-toc(clock)));
[pairs,report.screen]=ctocscreen.v3TimeCandidates(node,eph,screen,stream,pheromone,clock,coverage);
available=1:numel(pairs); selected=[];
if ~c.cost_guidance_enabled&&node.t==0&&isfield(node,'seed_target')&&node.seed_target>0
 matches=find(cellfun(@(p)p.target==node.seed_target,pairs));
 if ~isempty(matches)
  [~,k]=min(cellfun(@(p)abs(p.dt-node.seed_duration_s),pairs(matches)));
  selected=matches(k); available(available==selected)=[];
 end
end
if c.cost_guidance_enabled&&~isempty(available)
 [~,pick]=max(cellfun(@(p)p.weight,pairs)); selected=pick; available(available==pick)=[];
end
for k=1:min(c.guided_attempts-numel(selected),numel(available))
 choices=available;
 if c.cost_guidance_enabled&&~isempty(selected)
  targets=unique(cellfun(@(p)p.target,pairs(selected)));
  if numel(targets)<c.guided_distinct_targets
   other=choices(~cellfun(@(p)ismember(p.target,targets),pairs(choices)));
   if ~isempty(other), choices=other; end
  end
 end
 weights=cellfun(@(p)p.weight,pairs(choices));
 % Guard the selection: a degenerate or non-finite weight vector must not turn
 % pick into a non-scalar, which would invalidate the index downstream.
 if isempty(weights)||any(~isfinite(weights))||sum(weights)<=0
  pick=1;
 elseif rand(stream)<c.branch_greedy_probability
  [~,pick]=max(weights);
 else
  if c.cost_guidance_enabled, weights=weights/sum(weights);
  else, weights=(1-c.exploration)*weights/sum(weights)+c.exploration/numel(weights); end
  pick=find(cumsum(weights)>=rand(stream),1);
  if isempty(pick), pick=numel(weights); end
 end
 pick=pick(1);
 selected(end+1)=choices(pick); available(available==choices(pick))=[];
end
% Each selected pair receives a share before any local time polishing.
for k=1:numel(selected)
 if toc(clock)>=c.budget_s, break; end
 pair=pairs{selected(k)}; report.selected_pairs{end+1}=pair;
 allowance=.65*max(0,c.budget_s-toc(clock))/(numel(selected)-k+1);
 trial(pair.target,pair.dt,toc(clock)+allowance,'coarse');
end
if isempty(actions), return; end
[~,order]=sort(cellfun(@(a)norm(a.dv),actions)); seeds=actions(order);
for k=1:min(c.guided_attempts,numel(seeds))
 best=seeds{k}; step=max(c.action_min_duration_s,best.dt*c.time_refine_fraction);
 for it=1:c.time_refine_iterations
  for sign=[-1 1]
   if toc(clock)>=c.budget_s, break; end
   dt=min(eph.model.horizon_s-node.t,max(c.action_min_duration_s,best.dt+sign*step));
   before=numel(actions); trial(best.target,dt,c.budget_s,'local_time');
   if numel(actions)>before
    costs=cellfun(@(a)norm(a.dv),actions(before+1:end)); [cost,index]=min(costs);
    if cost<norm(best.dv), best=actions{before+index}; report.local_improvements=report.local_improvements+1; end
   end
  end
  step=step/2;
 end
end

 function trial(target,dt,limit,stage)
  if any(cellfun(@(r)r.target==target&&abs(r.dt-dt)<c.action_min_duration_s,report.trials)), return; end
  item=struct('target',target,'dt',dt,'stage',stage,'passed',false,'branches',0,'dv_km_s',Inf,'failure','', ...
   'branch_failures',{{}});
  try
   local=c; local.budget_s=min(c.budget_s,limit);
   goal=ctocscreen.v3QueryTargets(eph,target,node.t+dt).';
   departureSeed=[];
   if strcmp(stage,'coarse')&&isfield(pair,'v_depart'), departureSeed=pair.v_depart; end
   [~,~,branches]=ctocscreen.v3GuidedTransfer(node.state,node.t,node.t+dt,goal,eph.model,local,clock,departureSeed);
   for b=1:numel(branches)
    a=branches{b}; a.target=target; a.dt=dt; a.departure_s=node.t; a.stage=stage;
    feedback=ctocscreen.v3PolicyFeedback(pheromone,node,target,dt,a.dv,c);
    report.policy_feedback{end+1}=feedback;
    if feedback.total_dv_km_s>c.search_max_dv_km_s*c.rejected_proposal_factor
     item.failure='rejected_proposal_envelope_exceeded'; continue;
    end
    try
     % Charge visit scanning before further time-search trials, and reuse it
     % when packaging children. Finish the already integrated branch batch.
     hints=nan(35,1); hints(target)=node.t+dt;
     a.arc.scan=ctocscreen.v3ScanArc(a.arc.solution, ...
      @(ids,t,mode)ctocscreen.v3QueryTargets(eph,ids,t,mode),c,hints,find(~node.visited));
     actions{end+1}=a;
     item.branches=item.branches+1; item.dv_km_s=min(item.dv_km_s,norm(a.dv));
    catch err
     item.branch_failures{end+1}=struct('id',err.identifier,'message',err.message,'branch',a.branch);
    end
   end
   item.passed=item.branches>0;
  catch err
   item.failure=[err.identifier ': ' err.message];
  end
  report.trials{end+1}=item;
 end
end
