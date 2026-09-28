function [pairs,report]=v3CostCandidates(node,eph,c,stream,pheromone,clock)
%V3COSTCANDIDATES Plane-crossing/time proposals ranked by departure impulse.
% Lambert is a seed only. Acceptance always uses actual J2 propagation.
ids=find(~node.visited); pairs={};
report=struct('propagations',0,'pairs_evaluated',0,'failures',{{}}, ...
 'method','plane_windows_lambert_departure_cost','pheromone_queries',0, ...
 'nonneutral_queries',0,'prefix_guidance_hits',0,'negative_guidance_hits',0, ...
 'estimated_cost_rejections',0,'lambert_calls',0,'time_refinements',0,'natural_windows',0, ...
 'family_separation_hits',0,'family_completion_hits',0);
if isempty(ids), return; end
% Target families: visiting the same altitude band and orbital plane is far
% cheaper than hopping between families, so the ranking separates a leg's own
% cost from the family round-trip it commits to (Izzo et al. phasing indicators
% are used the same way: rank opportunities cheaply, verify with the real model).
[families,homeBand,homeInc,familyTargets,familyOf]=localfamilies(eph,node,ids);
stageStart=toc(clock); enumerationEnd=stageStart+c.cost_enumeration_fraction*max(0,c.budget_s-stageStart);
timeLeft=eph.model.horizon_s-node.t; slot=timeLeft/numel(ids);
look=min(timeLeft,c.cost_lookahead_s);
dtgrid=unique(min(look,[1800 3600 7200 c.cost_time_budget_factor*slot ...
 linspace(max(1800,.15*slot),look,c.cost_time_samples)]));
% Far-field phasing: cheap legs often need a long coast, but widening the whole
% lookahead makes every enumeration expensive. Instead keep the near field dense
% and add a few coarse far-field times, letting the task-time pressure term
% (v3LambertScore) reject the ones that would not leave time for the remaining
% targets. This is a proposal grid, not a time limit.
if c.far_field_samples>0
 farHorizon=min(timeLeft,max(look,c.far_field_horizon_s));
 farTimes=min(timeLeft,max(1800,farHorizon*(1:c.far_field_samples)/(c.far_field_samples+1)));
 dtgrid=unique([dtgrid(:);farTimes(:)]);
end
grid=linspace(node.t,node.t+look,max(3,ceil(look/1200)+1));
positions=ctocscreen.v3QueryTargets(eph,ids,grid,'grid');
natural=[];
try
 [~,~,natural]=ctocscreen.v3Arc(node.state,node.t,node.t+look,eph.model,c);
 naturalStates=deval(natural,grid);
 report.propagations=report.propagations+1;
catch
 % Failed lookahead is not a negative observation about a target.
end
normal=cross(node.state(1:3),node.state(4:6)); normal=normal/norm(normal);
allTimes=cell(numel(ids),1);
for k=1:numel(ids)
 if toc(clock)>=enumerationEnd, break; end
 r=reshape(positions(k,:,:),3,[]).'; z=r*normal(:);
 crossing=find(z(1:end-1).*z(2:end)<=0);
 ts=dtgrid; nearTimes=[];
 if ~isempty(natural)
  miss=vecnorm(r-naturalStates(1:3,:).',2,2);
  minima=find(miss(2:end-1)<=miss(1:end-2)&miss(2:end-1)<=miss(3:end))+1;
  [~,sorted]=sort(miss(minima));
  for jj=minima(sorted(1:min(c.cost_natural_windows,numel(sorted)))).'
   if toc(clock)>=enumerationEnd, break; end
   try
    t=fminbnd(@naturalMiss,grid(jj-1),grid(jj+1),optimset('Display','off','MaxFunEvals',15,'TolX',1));
    nearTimes(end+1)=t-node.t; report.natural_windows=report.natural_windows+1;
   catch
   end
  end
 end
 for j=crossing(:).'
  if toc(clock)>=enumerationEnd, break; end
  try
   t=fzero(@(t)dot(ctocscreen.v3QueryTargets(eph,ids(k),t),normal),grid(j:j+1));
   ts=[ts,t-node.t+[-600 0 600]];
  catch
   % A failed window root is not evidence against the target.
  end
 end
 if isempty(node.schedule.maneuver_times_s)&&isfield(node,'seed_target')&&node.seed_target==ids(k)
  ts=[ts node.seed_duration_s];
 end
 ts=unique(ts(ts>=c.action_min_duration_s&ts<=look));
 % Far-field times were sampled outside the near lookahead; keep them as long as
 % they stay inside the mission horizon.
 if c.far_field_samples>0
  ts=unique([ts(:);min(timeLeft,dtgrid(dtgrid>look))]);
  ts=ts(ts>=c.action_min_duration_s);
 end
 % Cover middle/long times before the short-time tail if the stage expires.
 indices=unique([ceil(numel(ts)/2),numel(ts),1,ceil(numel(ts)/4),ceil(3*numel(ts)/4),1:numel(ts)],'stable');
 indices=indices(indices>=1&indices<=numel(ts)); ts=ts(indices);
 if isempty(node.schedule.maneuver_times_s)&&isfield(node,'seed_target')&&node.seed_target==ids(k)
  [~,hint]=min(abs(ts-node.seed_duration_s)); ts=ts([hint setdiff(1:numel(ts),hint,'stable')]);
 end
 allTimes{k}=unique([nearTimes,ts],'stable');
 if isfield(node,'guidance_times_s')
  hint=node.guidance_times_s(ids(k))-node.t;
  if isfinite(hint)&&hint>=c.action_min_duration_s&&hint<=look
   allTimes{k}=unique([hint,allTimes{k}],'stable');
  end
 end
end
% Interleave targets so a wall deadline cannot always exclude the last IDs.
order=randperm(stream,numel(ids));
if isfield(node,'priority_targets')
 focus=ismember(ids(order),node.priority_targets);
 order=[order(focus),order(~focus)];
end
for round=1:max(cellfun(@numel,allTimes))
 for k=order
  if toc(clock)>=enumerationEnd, break; end
  if round>numel(allTimes{k}), continue; end
  dt=allTimes{k}(round);
  if toc(clock)>=enumerationEnd, break; end
  try
   [best,chosen,raw,plane]=ctocscreen.v3LambertScore(node,ids(k),dt,eph,c);
   report.lambert_calls=report.lambert_calls+1;
   if isempty(chosen), continue; end
   % Estimated cost never causes negative learning; only a corrected impulse does.
   if node.J+raw>c.search_max_dv_km_s*c.rejected_proposal_factor*1.15
    report.estimated_cost_rejections=report.estimated_cost_rejections+1; continue;
   end
   key=ctocscreen.v3StateKey(node,3,ids(k),dt,c);
   [exact,~,a]=ctocscreen.v3Pheromone('get',pheromone,key,0,c);
   policyKey=ctocscreen.v3PolicyKey(node,3,ids(k),dt);
   [backoff,~,b]=ctocscreen.v3Pheromone('get',pheromone,policyKey,0,c);
   ph=exact*backoff/c.pheromone_baseline;
   report.pheromone_queries=report.pheromone_queries+1;
   report.nonneutral_queries=report.nonneutral_queries+(abs(ph-c.pheromone_baseline)>1e-12);
   report.prefix_guidance_hits=report.prefix_guidance_hits+(a.prefix_increment+b.prefix_increment>0);
   report.negative_guidance_hits=report.negative_guidance_hits+(b.negative_penalty>0);
   pairs{end+1}=struct('target',ids(k),'dt',dt,'estimate',best,'estimated_dv_km_s',raw, ...
    'estimated_inclination_change_deg',plane.inclination_change_deg, ...
    'pheromone',ph,'natural_miss_km',NaN,'weight',0,'v_depart',chosen, ...
    'family',families(k),'family_cost_km_s',0);
  catch err
   report.failures{end+1}=struct('id',err.identifier,'message',err.message,'dt',dt);
  end
 end
 if toc(clock)>=enumerationEnd, break; end
end
report.pairs_evaluated=numel(pairs);
if isempty(pairs), return; end
scores=cellfun(@(p)p.estimate,pairs); [~,rank]=sort(scores);
refined=0; selectedTargets=[];
for index=rank
 if refined>=c.cost_refine_count||toc(clock)>=c.budget_s, break; end
 original=pairs{index};
 if ismember(original.target,selectedTargets), continue; end
 selectedTargets(end+1)=original.target; refined=refined+1;
 width=max(600,.15*original.dt); lo=max(c.action_min_duration_s,original.dt-width); hi=min(look,original.dt+width);
 try
  t=fminbnd(@objective,lo,hi,optimset('Display','off','MaxFunEvals',18,'TolX',1));
  [score,v,dv,p]=ctocscreen.v3LambertScore(node,original.target,t,eph,c);
  if score<original.estimate&&~isempty(v)
   report.time_refinements=report.time_refinements+1;
   item=original; item.dt=t; item.estimate=score; item.estimated_dv_km_s=dv;
   item.estimated_inclination_change_deg=p.inclination_change_deg; item.v_depart=v;
   key=ctocscreen.v3StateKey(node,3,item.target,t,c);
   exact=ctocscreen.v3Pheromone('get',pheromone,key,0,c);
   backoff=ctocscreen.v3Pheromone('get',pheromone,ctocscreen.v3PolicyKey(node,3,item.target,t),0,c);
   item.pheromone=exact*backoff/c.pheromone_baseline;
   % Store one struct per cell element: appending the struct array itself would
   % turn the cell into a 1-by-N struct and break every later field access.
   pairs=[pairs,num2cell(item).'];  %#ok<AGROW>
  end
 catch
  % Retain the coarse seed if local time refinement cannot improve it.
 end
end
scores=cellfun(@(p)p.estimate,pairs); minimum=min(scores);
scale=max(.05,min(.5,(c.search_max_dv_km_s-node.J)/max(1,numel(ids))));
for k=1:numel(pairs)
 boost=1;
 if isfield(node,'priority_targets')&&ismember(pairs{k}.target,node.priority_targets), boost=2; end
 targetId=pairs{k}.target;
 familyId=familyOf(targetId);
 memberIds=familyTargets{targetId};
 % Family separation: a leg that leaves the current altitude band or plane
 % commits the tour to a return trip. This is a ranking term only; the real
 % cost stays the propagated J2 impulse and the acceptance gate is unchanged.
 sep=ctocscreen.v3FamilySeparation(node.t,pairs{k}.dt,eph,pairs{k}.target,c,homeBand,homeInc);
 % Family completion: the last targets of a nearly finished family are the
 % cheapest way to stop paying that round trip; PHASING-INCREMENT rewards legs
 % cheaper than the running tour average (the cheap insertions that make the
 % phasing increments) instead of the expensive family transitions.
 familyReward=0;
 if numel(memberIds)>1&&familyId==familyOf(targetId)
  remainingInFamily=sum(~node.visited(memberIds));
  if remainingInFamily<=max(1,ceil(c.family_completion_window*numel(memberIds)))
   familyReward=c.family_completion_weight*(1+numel(memberIds)-remainingInFamily)/numel(memberIds);
  end
 end
 pursuit=max(0,median(scores)-pairs{k}.estimate)*c.phasing_pursuit_weight;
 pairs{k}.family_cost_km_s=sep;
 if sep>0, report.family_separation_hits=report.family_separation_hits+1; end
 if familyReward>0, report.family_completion_hits=report.family_completion_hits+1; end
 pairs{k}.weight=boost*max(0,1+c.family_penalty_weight*(familyReward-sep)+pursuit)* ...
  pairs{k}.pheromone^c.pheromone_alpha*exp(-min(50,(scores(k)-minimum)/scale));
end
end

function [family,homeBand,homeInc,byFamily,familyOf]=localfamilies(eph,node,ids)
%LOCALFAMILIES Altitude-band and plane groups of the targets still to visit.
% Search guidance only: it never enters the cost, the constraints or a result.
states=eph.states0(ids,:); r=vecnorm(states(:,1:3),2,2);
invA=2./r-sum(states(:,4:6).^2,2)/eph.model.mu; a=1./invA;
h=cross(states(:,1:3),states(:,4:6),2); inc=acosd(h(:,3)./vecnorm(h,2,2));
band=round(a/2.5e3); plane=round(inc/10);
[~,~,family]=unique([band(:) plane(:)],'rows');
homeBand=round(norm(node.state(1:3))/2.5e3);
hh=cross(node.state(1:3),node.state(4:6)); homeInc=acosd(hh(3)/max(norm(hh),eps));
byFamily=cell(35,1); familyOf=zeros(35,1);
for j=1:numel(ids)
 idx=find(family==family(j));
 byFamily{ids(j)}=ids(idx).';
 familyOf(ids(j))=family(j);
end
end
 function value=objective(t)
  value=1e6;
  if toc(clock)>=c.budget_s, return; end
  try, value=ctocscreen.v3LambertScore(node,original.target,t,eph,c); catch, end
  if ~isfinite(value), value=1e6; end
 end
 function value=naturalMiss(t)
  x=deval(natural,t); target=ctocscreen.v3QueryTargets(eph,ids(k),t);
  value=sum((x(1:3).'-target).^2);
 end
