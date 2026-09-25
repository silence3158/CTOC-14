function probe=v3ReplanProbe(s,eph,c,budget)
%V3REPLANPROBE Fixed-control propagation and witness probes for construction.
% Unconfirmed hints do not prove absence; this is never a verification gate.
clock=tic; s=ctocscreen.v3Normalize(s,eph.model);
probe=struct('schedule_key',ctocscreen.v3ScheduleKey(s),'status','budget_exhausted', ...
 'distance_km',inf(35,1),'planning_distance_km',inf(35,1), ...
 'witness_times_s',nan(35,1),'plan_times_s',nan(35,1), ...
 'trace',struct('times',[],'states',[],'arcs',{{}},'height_passed',[]), ...
 'elapsed_s',0,'failure_reason','','method','fixed_controls_witness_probe');
plan=s.witness_times_s;
if isfield(s,'visit_plan_times_s'), plan=s.visit_plan_times_s; end
missing=~isfinite(plan); fallback=linspace(0,s.duration_s,37).';
fallback=fallback(2:36); plan(missing)=fallback(missing); probe.plan_times_s=plan;
hints=[s.witness_times_s,plan];
knots=unique([0;s.maneuver_times_s;s.duration_s]);
x=ctocscreen.initialState(s.initial_q,eph.model.mu,eph.model.re).';
local=c; local.shooting_reltol=min(c.shooting_reltol,1e-12);
try
 for k=1:numel(knots)-1
  if toc(clock)>=budget, probe.elapsed_s=toc(clock); return; end
  burn=find(s.maneuver_times_s==knots(k),1);
  if ~isempty(burn), x(4:6)=x(4:6)+s.delta_v_km_s(burn,:).'; end
  probe.trace.times(end+1,1)=knots(k); probe.trace.states(end+1,:)=x.';
  [x,~,arc]=ctocscreen.v3Arc(x,knots(k),knots(k+1),eph.model,local);
  probe.trace.arcs{end+1}=arc;
  height=ctocscreen.v3Height(arc,eph.model,local);
  probe.trace.height_passed(end+1)=height.passed;
  for col=1:size(hints,2)
   ids=find(isfinite(hints(:,col))&hints(:,col)>=knots(k)&hints(:,col)<=knots(k+1));
   if isempty(ids), continue; end
   times=hints(ids,col); y=deval(arc,times);
   goals=ctocscreen.v3QueryTargets(eph,ids,times,'pairs');
   distance=vecnorm(y(1:3,:).'-goals,2,2); better=distance<probe.distance_km(ids);
   if col==2, probe.planning_distance_km(ids)=distance; end
   probe.distance_km(ids(better))=distance(better);
   probe.witness_times_s(ids(better))=times(better);
  end
 end
 probe.trace.times(end+1,1)=s.duration_s;
 burn=find(s.maneuver_times_s==s.duration_s,1);
 if ~isempty(burn), x(4:6)=x(4:6)+s.delta_v_km_s(burn,:).'; end
 probe.trace.states(end+1,:)=x.'; probe.status='complete';
catch err
 probe.status='propagation_failure'; probe.failure_reason=[err.identifier ': ' err.message];
end
probe.elapsed_s=toc(clock);
end
