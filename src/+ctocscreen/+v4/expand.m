function [children,report,memory]=expand(parent,eph,c,stream,memory,budget)
%EXPAND Coasts and delayed impulse proposals from the current real trajectory.
clock=tic; children={}; m=eph.model; q=parent.q; t0=q.T; x=parent.actual.final_state;
remaining=find(parent.actual.distance_km>1); report=struct('enumerated',0,'guided',0, ...
 'children',0,'seconds',0,'coast_children',0,'failed_guidance',0,'errors',{{}});
if isempty(remaining)||t0>=m.horizon_s, return; end
available=m.horizon_s-t0;
look=min([available,c.lookahead_s,3*available/max(1,numel(remaining))]);
if look<=1, return; end
% Natural continuation is a first-class child, even without a new visit.
if parent.zero_gain<c.max_zero_gain
 coast=q; coast.T=t0+look;
 [a,tr]=ctocscreen.v4.replay(coast,eph,c,false,parent);
 if a.initial_passed&&a.height_passed&&~strcmp(a.status,'propagation_failure')
  node=makeChild(parent,coast,a,tr,'coast'); children{end+1}=node; report.coast_children=1;
 end
end
if toc(clock)>=budget, report.seconds=toc(clock); report.children=numel(children); return; end
dt=unique([max(600,look*[.06 .12 .25 .45 .7 1]),min(look,parent.seed_duration)]);
dt=dt(dt>0&dt<=available); if numel(dt)>c.time_samples, dt=dt(round(linspace(1,numel(dt),c.time_samples))); end
% Actual target geometry only; keep exploration targets outside the nearest group.
h=cross(x(1:3),x(4:6)); h=h/max(norm(h),eps);
[rr,~]=ctocscreen.v3QueryTargets(eph,remaining,t0+min(look,median(dt)));
rankCost=abs(rr*h)+.05*abs(vecnorm(rr,2,2)-norm(x(1:3)));
[~,rank]=sort(rankCost); chosen=remaining(rank(1:min(c.target_count,numel(rank))));
if numel(remaining)>numel(chosen), chosen(end)=remaining(randi(stream,numel(remaining))); end
if t0==0&&~ismember(parent.seed_target,chosen), chosen(1)=parent.seed_target; end
chosen=unique(chosen,'stable');
pairs=zeros(numel(chosen)*numel(dt),2); at=0;
for it=1:numel(dt)
 for j=1:numel(chosen), at=at+1; pairs(at,:)=[j it]; end
end
pairs=pairs(randperm(stream,size(pairs,1)),:);
seeds=struct('target',{},'arrival',{},'departure',{},'x',{},'v',{},'score',{});
enumerationEnd=min(budget*.45,toc(clock)+max(.2,budget*.3));
for index=1:size(pairs,1)
 if toc(clock)>=enumerationEnd, break; end
 id=chosen(pairs(index,1)); duration=dt(pairs(index,2)); wait=0;
 if mod(index,4)==0, wait=min(duration*.25,3600); end
 if wait>0
  % Reintegrate from the last real impulse; a temporary coast endpoint is not a control event.
  if isempty(q.tau), anchor=0; xa=q.x0;
  else, anchor=q.tau(end); xa=ctocscreen.v4.stateAt(parent.trace,anchor,'post'); end
  xp=ctocscreen.v3Arc(xa,anchor,t0+wait,m,c,false,true);
 else, xp=x; end
 rt=ctocscreen.v3QueryTargets(eph,id,t0+duration);
 policy=struct('max_revolutions',c.max_revolutions,'endpoint_tol_km',.005);
 try
  branches=ctocscreen.v3LambertBranches(xp(1:3),rt,duration-wait,m.mu,policy);
  for b=1:numel(branches)
   dv=branches(b).v_depart(:)-xp(4:6);
   h1=cross(xp(1:3),xp(4:6)); h2=cross(xp(1:3),branches(b).v_depart(:));
   inc=abs(acosd(max(-1,min(1,h1(3)/norm(h1))))-acosd(max(-1,min(1,h2(3)/norm(h2)))));
   penalty=(max(0,inc-c.plane_threshold_deg)/c.plane_threshold_deg)^2;
   score=norm(dv)+c.plane_weight*min(4,penalty);
   seeds(end+1)=struct('target',id,'arrival',t0+duration,'departure',t0+wait, ...
    'x',xp,'v',branches(b).v_depart(:),'score',score); %#ok<AGROW>
  end
  report.enumerated=report.enumerated+1;
 catch err
  report.errors{end+1}=err.identifier;
 end
end
if ~isempty(seeds)
 [~,order]=sort([seeds.score]);
 % Preserve the cheapest proposal and sample an alternative from distinct windows.
 if numel(order)>2&&rand(stream)<c.exploration
  pick=randi(stream,[2 min(numel(order),10)]); order([2 pick])=order([pick 2]);
 end
 used=[];
 for k=order
  if toc(clock)>=budget||report.guided>=c.branch_count, break; end
  seed=seeds(k); signature=[seed.target,seed.arrival];
  if ~isempty(used)&&ismember(signature,used,'rows'), continue; end
  used(end+1,:)=signature; %#ok<AGROW>
  report.guided=report.guided+1;
  gc=ctocscreen.v3Defaults(struct('budget_s',max(.02,budget-toc(clock)), ...
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
    node=makeChild(parent,child,a,tr,'delayed_impulse');
    [memory,~]=ctocscreen.v4.feedback(memory,'observe',node,c);
    children{end+1}=node;
   end
  catch err
   report.failed_guidance=report.failed_guidance+1; report.errors{end+1}=err.identifier;
  end
 end
end
report.children=numel(children); report.seconds=toc(clock);
end
function child=makeChild(parent,q,a,tr,origin)
child=parent; child.q=tr.q; child.q.witness=a.witness_times_s;
child.actual=a; child.trace=tr; child.origin=origin; child.generation=parent.generation+1;
child.attempts=0; child.zero_gain=0;
if a.visit_count<=parent.actual.visit_count, child.zero_gain=parent.zero_gain+1; end
end
