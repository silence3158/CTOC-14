function [s,r]=rebuildSchedule(plan,p,cfg)
%REBUILDSCHEDULE Retarget suffix from actual states; coast through skip slots.
s=struct('schema_version','fragment_schedule_v2','dynamics_id','two_body', ...
 'initial_q',plan.initial_q,'maneuver_times_s',zeros(0,1),'delta_v_km_s',zeros(0,3), ...
 'duration_s',plan.times(end),'event_target_ids',plan.ids(:),'event_times_s',plan.times(:));
r=struct('passed',false,'total_dv_km_s',Inf,'status','infeasible');
try
 assert(all(diff([0;plan.times(:)])>0)&&plan.times(end)<=p.horizon_s);
 state=ctocscreen.initialState(plan.initial_q,p.mu_km3_s2); t=0;
 for k=1:numel(plan.ids)
  dt=plan.times(k)-t; w=plan.waits(k); assert(w>=0&&w<dt);
  a=ctocscreen.checkArc(state,w,p); assert(strcmp(a.status,'ok')&&a.min_altitude_km>=200);
  state=ctocscreen.propagateTwoBody(state,w,p.mu_km3_s2);
  target=ctocscreen.targetStates(p,plan.ids(k),plan.times(k));
  if ~plan.skip(k)
   if all(isfinite(plan.fixed_dv(k,:)))
    u=plan.fixed_dv(k,:);
   else
    [bs,~]=ctocscreen.enumerateBranches(state(1:3),target(1:3),dt-w,p.mu_km3_s2,cfg.branch_policy);
    assert(~isempty(bs));
    if all(isfinite(plan.reference_v(k,:)))
     costs=arrayfun(@(b)norm(b.v_depart-plan.reference_v(k,:)),bs);
    else
     costs=arrayfun(@(b)norm(b.v_depart-state(4:6)),bs);
    end
    [~,ix]=sort(costs); found=false;
    for b=ix(:)'
     arc=ctocscreen.checkArc([state(1:3) bs(b).v_depart],dt-w,p);
     if strcmp(arc.status,'ok')&&arc.min_altitude_km>=200, found=true; break; end
    end
    assert(found); u=bs(b).v_depart-state(4:6);
   end
   s.maneuver_times_s(end+1,1)=t+w; s.delta_v_km_s(end+1,:)=u;
   state(4:6)=state(4:6)+u;
  end
  a=ctocscreen.checkArc(state,dt-w,p); assert(strcmp(a.status,'ok')&&a.min_altitude_km>=200);
  state=ctocscreen.propagateTwoBody(state,dt-w,p.mu_km3_s2);
  assert(norm(state(1:3)-target(1:3))<=cfg.fragment_tolerance_km);
  t=plan.times(k);
 end
 r=ctocscreen.propagateSchedule(s,p,false);
 r.build_completed=true; % Numeric trajectory exists even if initial-orbit bounds reject it.
catch err
 r.failure_reason=err.message;
end
end
