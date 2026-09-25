function [proposals,report]=v3Replan(original,eph,c,stream,pheromone,budget,attemptNumber,cutRound)
%V3REPLAN TAS/SRS-inspired free-flyby rebuilding, followed by full-task repair.
if nargin<7, attemptNumber=1; end
if nargin<8, cutRound=attemptNumber; end
clock=tic; proposals={};
rootKey=sprintf('%.17g/',original.initial_q);
if isfield(original,'root_key'), rootKey=original.root_key; end
report=struct('status','not_constructed','reason','','input_key',ctocscreen.v3ScheduleKey(original), ...
 'attempt_number',attemptNumber,'cut_round',cutRound,'root_key',rootKey, ...
 'focus_targets',[],'cut_times_s',[],'available_cut_times_s',[], ...
 'attempts',{{}},'elapsed_s',0,'budget_overrun_s',0,'replay_seconds',0, ...
 'prefix_seconds',0,'full_feasibility_verified',false,'target_ids',1:35);
try
 original=ctocscreen.v3Normalize(original,eph.model);
 costDriven=sum(vecnorm(original.delta_v_km_s,2,2))>c.search_max_dv_km_s;
 report.cost_driven=costDriven;
 probe=ctocscreen.v3ReplanProbe(original,eph,c,max(0,budget-toc(clock)));
 report.replay_seconds=probe.elapsed_s; report.probe_status=probe.status;
 if strcmp(probe.status,'budget_exhausted')
  report.reason=probe.status; report.probe_failure=probe.failure_reason;
 else
  missing=find(probe.distance_km>1|probe.planning_distance_km>1); report.focus_targets=missing.';
  report.unconfirmed_targets=find(probe.distance_km>1).';
  report.plan_miss_targets=find(probe.planning_distance_km>1).';
  if strcmp(probe.status,'propagation_failure')
   cuts=0; modes={'propagation_replan'}; report.probe_failure=probe.failure_reason;
  else
   [cuts,modes]=locations(original,probe,missing,eph.model,c);
  end
  report.available_cut_times_s=cuts;
  count=min(c.replan_cut_count,numel(cuts));
  % Rotate the first slot too, so an overrun cannot starve later entries.
  selected=mod((cutRound-1)+(0:count-1),numel(cuts))+1;
  report.cut_times_s=cuts(selected);
  for k=1:count
   if toc(clock)>=budget, report.reason='budget_exhausted'; break; end
   cut=cuts(selected(k)); mode=modes{selected(k)};
   a=struct('mode',mode,'cut_s',cut,'status','failed','reason','', ...
    'prefix_state',[],'released_target_ids',[],'completion',[], ...
    'failures',{{}},'proposals',0,'elapsed_s',0);
   attemptClock=tic;
   try
    prefixClock=tic; start=ctocscreen.v3PrefixNode(original,cut,eph,c,probe);
    report.prefix_seconds=report.prefix_seconds+toc(prefixClock);
    a.prefix_state=start.state; a.released_target_ids=find(~start.visited).';
    a.retained_dv_km_s=start.J;
    start.priority_targets=missing;
    if isempty(missing)
     later=probe.plan_times_s>cut;
     start.priority_targets=find(later&probe.plan_times_s<=nextBoundary(original,cut));
    end
    start.guidance_times_s=probe.plan_times_s;
    if all(start.visited)
     % A complete pre-burn prefix can remove a redundant terminal suffix.
     nodes={start};
    elseif toc(clock)<budget
     local=c; local.completion_seconds=max(.001,(budget-toc(clock))/(count-k+1));
     [~,a.completion,a.failures,nodes]=ctocscreen.v3Complete(start,eph,local,stream,pheromone,false, ...
      c.cold_completion_enabled&&(costDriven||sum(isfinite(original.witness_times_s))<35));
    else
     nodes={}; a.reason='budget_exhausted';
    end
    for j=1:numel(nodes)
     n=nodes{j};
     if n.t<cut||n.t<=0, continue; end
     if n.t==cut&&~all(n.visited), continue; end
     if costDriven
      % A new partial tail must be completed from its real state. Do not
      % append the rejected parent's controls or resurrect its witnesses.
      s=n.schedule; s.validation_level='proposal';
      if isfield(s,'visit_plan_times_s'), s=rmfield(s,'visit_plan_times_s'); end
     else
      s=compose(original,n,probe.plan_times_s,eph.model);
     end
     difference=ctocscreen.v3MutationDifference(original,s,c,eph.model);
     if ~difference.novel, continue; end
     key=ctocscreen.v3ScheduleKey(s);
     if any(cellfun(@(p)strcmp(ctocscreen.v3ScheduleKey(p.schedule),key),proposals)), continue; end
     s.replan_origin=struct('input_key',report.input_key,'attempt_number',attemptNumber, ...
      'mode',mode,'cut_s',cut,'rebuilt_end_s',n.t);
     proposals{end+1}=struct('schedule',s,'target_ids',1:35,'mode',mode, ...
      'cut_s',cut,'difference',difference,'constructed_visits',sum(n.visited));
     if costDriven, proposals{end}.node=n; end
     a.proposals=a.proposals+1;
    end
    a.status='attempted';
   catch err
    a.reason=[err.identifier ': ' err.message];
   end
   a.elapsed_s=toc(attemptClock); report.attempts{end+1}=a;
  end
 end
 if ~isempty(proposals), report.status='constructed_proposals'; end
 if toc(clock)>=budget&&isempty(report.reason), report.reason='budget_exhausted'; end
catch err
 report.reason=[err.identifier ': ' err.message];
end
report.elapsed_s=toc(clock); report.budget_overrun_s=max(0,report.elapsed_s-budget);
end

function [cuts,modes]=locations(s,probe,missing,m,c)
cuts=[]; modes={};
% Move upstream of the first unsafe arc before considering downstream cuts.
bad=find(~probe.trace.height_passed,1);
if ~isempty(bad)
 cuts=probe.trace.times(bad); modes={'height_replan'};
end
if ~isempty(missing)
 [~,order]=sort(max(probe.distance_km(missing),probe.planning_distance_km(missing)),'descend');
 for id=missing(order).'
  event=probe.plan_times_s(id); earlier=s.maneuver_times_s(s.maneuver_times_s<event);
  cut=0; if ~isempty(earlier), cut=earlier(max(1,numel(earlier)-1)); end
  cuts(end+1)=cut; modes{end+1}='missing_target_replan';
 end
 if s.duration_s<m.horizon_s-c.action_min_duration_s
  cuts(end+1)=s.duration_s; modes{end+1}='tail_addition';
 end
else
 [~,order]=sort(vecnorm(s.delta_v_km_s,2,2),'descend');
 cuts=[cuts,s.maneuver_times_s(order).'];
 modes=[modes,repmat({'expensive_connection_replan'},1,numel(order))];
end
cuts(end+1)=0; modes{end+1}='broad_replan';
[cuts,index]=unique(cuts,'stable'); modes=modes(index);
if sum(vecnorm(s.delta_v_km_s,2,2))>c.search_max_dv_km_s
 [costCuts,~]=ctocscreen.v3CostReplanCuts(s,c);
 spent=arrayfun(@(t)sum(vecnorm(s.delta_v_km_s(s.maneuver_times_s<t,:),2,2)),cuts);
 safe=spent<c.search_max_dv_km_s; cuts=cuts(safe); modes=modes(safe);
 cuts=[costCuts cuts]; modes=[repmat({'cost_prefix_replan'},size(costCuts)) modes];
 [cuts,index]=unique(cuts,'stable'); modes=modes(index);
end
end

function t=nextBoundary(s,cut)
later=s.maneuver_times_s(s.maneuver_times_s>cut);
if isempty(later), t=s.duration_s; else, t=later(min(2,numel(later))); end
end

function s=compose(original,node,plan,m)
s=node.schedule;
if all(node.visited)
 s.duration_s=node.t;
else
 suffix=original.maneuver_times_s>node.t;
 s.maneuver_times_s=[s.maneuver_times_s;original.maneuver_times_s(suffix)];
 s.delta_v_km_s=[s.delta_v_km_s;original.delta_v_km_s(suffix,:)];
 s.duration_s=max(original.duration_s,node.t);
end
plan=min(s.duration_s,max(0,plan));
witnessed=isfinite(s.witness_times_s); plan(witnessed)=s.witness_times_s(witnessed);
s.visit_plan_times_s=plan; s.validation_level='proposal';
s=ctocscreen.v3Normalize(s,m);
end
