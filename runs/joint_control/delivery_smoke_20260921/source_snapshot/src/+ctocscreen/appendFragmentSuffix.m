function [s,r]=appendFragmentSuffix(base,entryIndex,f,p,cfg)
%APPENDFRAGMENTSUFFIX Preserve prefix impulses, splice fragment, retarget suffix.
s=base; r=struct('passed',false,'total_dv_km_s',Inf);
try
 tin=0; if entryIndex>0, tin=base.event_times_s(entryIndex); end
 assert(f.passed&&abs(f.t_in_s-tin)<1e-8);
 keep=base.maneuver_times_s<tin;
 s.maneuver_times_s=[base.maneuver_times_s(keep);f.maneuver_times_s];
 s.delta_v_km_s=[base.delta_v_km_s(keep,:);f.delta_v_km_s];
 ids=base.event_target_ids; remaining=(entryIndex+1):numel(ids);
 remaining=remaining(~ismember(ids(remaining),f.event_target_ids));
 s.event_target_ids=[ids(1:entryIndex);f.event_target_ids;ids(remaining)];
 s.event_times_s=[base.event_times_s(1:entryIndex);f.event_times_s; base.event_times_s(remaining)];
 assert(all(diff(s.event_times_s)>0));
 old=ctocscreen.propagateSchedule(base,p,false);assert(old.passed);
 x=f.state_out; time=f.t_out_s;
 for j=remaining
  arrival=base.event_times_s(j); dt=arrival-time; assert(dt>60);
  oldStart=0;if j>1,oldStart=base.event_times_s(j-1);end
  burn=find(base.maneuver_times_s>=oldStart-1e-8 & base.maneuver_times_s<arrival,1,'last');
  wait=0;ref=[];
  if ~isempty(burn)
   wait=min(max(0,base.maneuver_times_s(burn)-oldStart),dt-60);
   ref=old.preburn_states(burn,4:6)+base.delta_v_km_s(burn,:);
  end
  coast=ctocscreen.checkArc(x,wait,p);assert(strcmp(coast.status,'ok')&&coast.min_altitude_km>=200);
  x=ctocscreen.propagateTwoBody(x,wait,p.mu_km3_s2);time=time+wait;dt=arrival-time;
  target=ctocscreen.targetStates(p,ids(j),arrival);
  bs=ctocscreen.enumerateBranches(x(1:3),target(1:3),dt,p.mu_km3_s2,cfg.branch_policy);
  assert(~isempty(bs)); if isempty(ref),ref=x(4:6);end
  [~,ii]=sort(arrayfun(@(b)norm(b.v_depart-ref),bs)); found=false;
  for b=ii(:)'
   a=ctocscreen.checkArc([x(1:3) bs(b).v_depart],dt,p);
   if strcmp(a.status,'ok')&&a.min_altitude_km>=200, found=true;break;end
  end
  assert(found); u=bs(b).v_depart-x(4:6);
  s.maneuver_times_s(end+1,1)=time; s.delta_v_km_s(end+1,:)=u;
  x(4:6)=x(4:6)+u; x=ctocscreen.propagateTwoBody(x,dt,p.mu_km3_s2); time=arrival;
 end
 s.duration_s=max(s.event_times_s); r=ctocscreen.propagateSchedule(s,p,false);
catch err
 r.failure_reason=err.message;
end
end
