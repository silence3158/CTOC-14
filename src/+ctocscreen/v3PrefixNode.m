function node=v3PrefixNode(schedule,cut,eph,c,probe)
%V3PREFIXNODE Extract a physically replayed prefix before the burn at cut.
assert(cut>=0&&cut<=schedule.duration_s);
s=schedule; if isfield(s,'visit_plan_times_s'), s=rmfield(s,'visit_plan_times_s'); end
keep=s.maneuver_times_s<cut; s.maneuver_times_s=s.maneuver_times_s(keep); s.delta_v_km_s=s.delta_v_km_s(keep,:);
s.witness_times_s(:)=NaN;
if cut==0
 s.duration_s=min(1,eph.model.horizon_s); x=ctocscreen.initialState(s.initial_q,eph.model.mu,eph.model.re);
 visited=false(35,1); J=0; penalty=0;
elseif nargin>=5
 assert(strcmp(probe.schedule_key,ctocscreen.v3ScheduleKey(schedule))&&strcmp(probe.status,'complete'), ...
  'ctocscreen:v3:replanProbe','Probe does not describe this trajectory.');
 s.duration_s=cut; arc=find(probe.trace.times<cut,1,'last');
 assert(~isempty(arc)&&all(probe.trace.height_passed(1:arc)), ...
  'ctocscreen:v3:replanPrefix','Unsafe prefix.');
 y=deval(probe.trace.arcs{arc},cut); x=y(1:6).';
 visited=probe.distance_km<=1 & probe.witness_times_s<=cut;
 s.witness_times_s(visited)=probe.witness_times_s(visited);
 J=sum(vecnorm(s.delta_v_km_s,2,2));
 penalty=ctocscreen.v3TracePlanePenalty(s,probe.trace,c);
else
 s.duration_s=cut; [r,trace]=ctocscreen.v3Replay(s,eph,c,false);
 assert(~strcmp(r.status,'propagation_failure')&&r.height_passed,'ctocscreen:v3:replanPrefix','Unsafe prefix.');
 x=r.final_state; visited=r.distance_km<=1; s.witness_times_s(visited)=r.witness_times_s(visited); J=r.total_dv_km_s;
 penalty=ctocscreen.v3TracePlanePenalty(s,trace,c);
end
if ~isfield(s,'root_key'), s.root_key=sprintf('%.17g/',s.initial_q); end
node=struct('schedule',s,'t',cut,'state',x,'visited',visited,'J',J,'keys',{{}}, ...
 'last_gain',0,'estimate',J,'seed_target',0,'seed_duration_s',NaN,'root_method','directed_replan');
node.inclination_penalty=penalty; node.estimate=J+c.plane_penalty_km_s*penalty;
end
