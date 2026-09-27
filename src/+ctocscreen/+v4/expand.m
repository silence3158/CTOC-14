function [children,report,memory]=expand(parent,eph,c,stream,memory,budget,options)
%EXPAND Children of the current real trajectory (A layer).
% Main source: timed Lambert transfers to target plane-crossing events
% (events.m, crossing.m). Old any-time grid proposals keep a minority share
% for diversity. All children are J2-corrected and replayed with fixed controls.
% Optional options (defaults keep the original behavior):
%  max_leg_dv       proposals whose two-body leg cost exceeds it are dropped
%                   before J2 guidance (branch-and-bound successor test with a
%                   predicted edge cost; the replayed child is still re-checked)
%  exclude          rows [target arrival_s] already tried at this node; same
%                   target within exclude_window_s counts as the same encounter
if nargin<7||isempty(options), options=struct(); end
maxLeg=Inf; if isfield(options,'max_leg_dv'), maxLeg=options.max_leg_dv; end
tried=zeros(0,2); if isfield(options,'exclude'), tried=options.exclude; end
window=600; if isfield(options,'exclude_window_s'), window=options.exclude_window_s; end
clock=tic; children={}; m=eph.model; q=parent.q; t0=q.T;
baseBudget=budget*(1-c.shared_enabled*c.shared_fraction);
remaining=find(parent.actual.distance_km>1); report=struct('enumerated',0,'guided',0, ...
 'children',0,'seconds',0,'coast_children',0,'failed_guidance',0,'errors',{{}}, ...
 'crossing_events',0,'crossing_proposals',0,'grid_proposals',0,'heuristic_H',NaN,'crossing_children',0, ...
 'budget_filtered',0,'excluded',0,'tried',zeros(0,2),'shared_calls',0,'shared_children',0, ...
 'shared_seconds',0,'shared_reports',{{}},'repair_seed_count',0);
if isempty(remaining)||t0>=m.horizon_s, return; end
[ev,info]=ctocscreen.v4.events(parent,eph,c);
report.crossing_events=numel(ev); report.heuristic_H=info.H;
seeds=ctocscreen.v4.crossing(parent,eph,c,ev,info,stream,.4*budget);
report.crossing_proposals=numel(seeds);
if rand(stream)<c.grid_share&&toc(clock)<.6*budget
 grid=gridSeeds(parent,eph,c,stream,remaining,.6*budget-toc(clock));
 report.grid_proposals=numel(grid); report.enumerated=numel(grid);
 seeds=[seeds,grid];
end
if ~isempty(seeds)
 keep=[seeds.dv]<=maxLeg;
 % One near-over-bound seed may still be improved by joint correction.
 if c.shared_enabled&&isfinite(maxLeg)
  repair=find(~keep&[seeds.dv]<=maxLeg+c.shared_seed_margin);
  if ~isempty(repair), [~,rp]=min([seeds(repair).dv]); keep(repair(rp))=true; report.repair_seed_count=1; end
 end
 report.budget_filtered=sum(~keep);
 for k=find(keep)
  if any(tried(:,1)==seeds(k).target&abs(tried(:,2)-seeds(k).arrival)<window)
   keep(k)=false; report.excluded=report.excluded+1;
  end
 end
 seeds=seeds(keep);
end
if ~isempty(seeds)
 % Sampling after spec eq. (30) with neutral pheromone: weight ~ eta^beta,
 % eta = 1/score, plus a uniform epsilon share (exploration).
 score=[seeds.score]; w=(1./max(score,1e-3)).^c.heuristic_beta; w=w/sum(w);
 prob=(1-c.exploration)*w+c.exploration/numel(w);
 used=[]; tries=0;
 while report.guided<c.branch_count&&tries<4*c.branch_count&&toc(clock)<baseBudget&&any(prob>0)
  tries=tries+1;
  % The cheapest proposal is always tried first; later picks are sampled.
  if report.guided==0, [~,k]=min(score); else, k=find(rand(stream)<=cumsum(prob)/sum(prob),1); end
  prob(k)=0; seed=seeds(k); signature=[seed.target,round(seed.arrival)];
  if ~isempty(used)&&ismember(signature,used,'rows'), continue; end
  used(end+1,:)=signature; %#ok<AGROW>
  report.guided=report.guided+1; report.tried(end+1,:)=[seed.target,seed.arrival];
  gc=ctocscreen.v3Defaults(struct('budget_s',max(.02,baseBudget-toc(clock)), ...
   'guided_branches',1,'guided_max_revolutions',c.max_revolutions,'search_radius_km',c.search_radius_km));
  gc.shooting_reltol=c.shooting_reltol; gc.verify_reltol=c.verify_reltol; gc.max_step_s=c.max_step_s;
  gc.height_margin_km=c.height_margin_km;
  try
   gcClock=tic; rt=ctocscreen.v3QueryTargets(eph,seed.target,seed.arrival);
   dv=ctocscreen.v3GuidedTransfer(seed.x,seed.departure,seed.arrival,rt,m,gc,gcClock,seed.v);
   child=q; child.tau(end+1,1)=seed.departure; child.u(end+1,:)=dv.'; child.T=seed.arrival;
   child.witness(seed.target)=seed.arrival;
   [a,tr]=ctocscreen.v4.replay(child,eph,c,false,parent);
   if a.initial_passed&&a.height_passed&&~strcmp(a.status,'propagation_failure')
    origin='delayed_impulse'; if strcmp(seed.kind,'crossing'), origin='crossing_transfer'; end
    node=makeChild(parent,child,a,tr,origin);
    [memory,~]=ctocscreen.v4.feedback(memory,'observe',node,c);
    children{end+1}=node; %#ok<AGROW>
    if strcmp(seed.kind,'crossing'), report.crossing_children=report.crossing_children+1; end
   end
  catch err
   report.failed_guidance=report.failed_guidance+1; report.errors{end+1}=err.identifier;
  end
 end
end
% Shared-arc proposals precede the caller's true-cost bound and beam pruning.
% Keep the original single-target child; both are alternatives of parent.
transferCount=numel(children);
for ks=1:transferCount
 if ~c.shared_enabled||toc(clock)>=budget-.05, break; end
 share=min(c.absorb_seconds,(budget-toc(clock))/(transferCount-ks+1));
 bt=tic; [paired,sr]=ctocscreen.v4.sharedArc(children{ks},parent,eph,c,share);
 report.shared_seconds=report.shared_seconds+toc(bt); report.shared_calls=report.shared_calls+1;
 report.shared_children=report.shared_children+numel(paired); report.shared_reports{end+1}=sr;
 for kp=1:numel(paired), [memory,~]=ctocscreen.v4.feedback(memory,'observe',paired{kp},c); end
 children=[children,paired]; %#ok<AGROW>
end
% A natural continuation stays available, but only when no transfer child exists:
% waiting is already represented by delayed departures along the coast.
if isempty(children)&&parent.zero_gain<c.max_zero_gain&&toc(clock)<budget
 look=min([m.horizon_s-t0,c.lookahead_s]);
 if look>1
  coast=q; coast.T=t0+look;
  [a,tr]=ctocscreen.v4.replay(coast,eph,c,false,parent);
  if a.initial_passed&&a.height_passed&&~strcmp(a.status,'propagation_failure')
   children{end+1}=makeChild(parent,coast,a,tr,'coast'); report.coast_children=1;
  end
 end
end
report.children=numel(children); report.seconds=toc(clock);
end
function seeds=gridSeeds(parent,eph,c,stream,remaining,budget)
% Pre-2026-09-26 any-time target/duration grid, kept as a diversity source.
clock=tic; m=eph.model; q=parent.q; t0=q.T; x=parent.actual.final_state;
seeds=struct('target',{},'departure',{},'arrival',{},'x',{},'v',{},'dv',{},'score',{},'kind',{});
available=m.horizon_s-t0;
look=min([available,c.lookahead_s,3*available/max(1,numel(remaining))]);
if look<=1, return; end
dt=unique([max(600,look*[.06 .12 .25 .45 .7 1]),min(look,parent.seed_duration)]);
dt=dt(dt>0&dt<=available); if numel(dt)>c.time_samples, dt=dt(round(linspace(1,numel(dt),c.time_samples))); end
h=cross(x(1:3),x(4:6)); h=h/max(norm(h),eps);
rr=ctocscreen.v3QueryTargets(eph,remaining,t0+min(look,median(dt)));
[~,rank]=sort(abs(rr*h)+.05*abs(vecnorm(rr,2,2)-norm(x(1:3))));
chosen=remaining(rank(1:min(c.target_count,numel(rank))));
if numel(remaining)>numel(chosen), chosen(end)=remaining(randi(stream,numel(remaining))); end
if t0==0&&~ismember(parent.seed_target,chosen), chosen(1)=parent.seed_target; end
chosen=unique(chosen,'stable');
policy=struct('max_revolutions',c.max_revolutions,'endpoint_tol_km',.005);
for id=chosen(:).'
 for d=dt
  if toc(clock)>=budget, return; end
  rt=ctocscreen.v3QueryTargets(eph,id,t0+d);
  try, b=ctocscreen.v3LambertBranches(x(1:3),rt,d,m.mu,policy); catch, continue; end
  for j=1:numel(b)
   v=b(j).v_depart(:); h1=cross(x(1:3),x(4:6)); h2=cross(x(1:3),v);
   inc=abs(acosd(max(-1,min(1,h1(3)/norm(h1))))-acosd(max(-1,min(1,h2(3)/norm(h2)))));
   penalty=(max(0,inc-c.plane_threshold_deg)/c.plane_threshold_deg)^2; dv=norm(v-x(4:6));
   seeds(end+1)=struct('target',id,'departure',t0,'arrival',t0+d,'x',x,'v',v, ...
    'dv',dv,'score',dv+c.plane_weight*min(4,penalty),'kind','grid'); %#ok<AGROW>
  end
 end
end
end
function child=makeChild(parent,q,a,tr,origin)
child=parent; child.q=tr.q; child.q.witness=a.witness_times_s;
child.actual=a; child.trace=tr; child.origin=origin; child.generation=parent.generation+1;
child.attempts=0; child.zero_gain=0; child.heuristic_H=NaN;
if a.visit_count<=parent.actual.visit_count, child.zero_gain=parent.zero_gain+1; end
end
