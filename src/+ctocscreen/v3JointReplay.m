function [out,trace]=v3JointReplay(s,eph,c,independent)
%V3JOINTREPLAY Keep a proposed visit assignment distinct from detected visits.
[out,trace]=ctocscreen.v3Replay(s,eph,c,independent);
out.physical_passed=out.passed;
if ~isfield(s,'visit_plan_times_s')||strcmp(out.status,'propagation_failure'), return; end
out.planning_distance_km=inf(35,1);
for id=c.joint_target_ids(:).'
 t=s.visit_plan_times_s(id);
 if ~isfinite(t), continue; end
 knot=find(trace.times==t,1);
 if isempty(knot)
  arc=find(trace.times<t,1,'last'); y=deval(trace.arcs{arc},t); x=y(1:3).';
 else
  x=trace.states(knot,1:3);
 end
 r=ctocscreen.v3QueryTargets(eph,id,t);
 out.planning_distance_km(id)=norm(x-r);
end
if any(out.planning_distance_km(c.joint_target_ids)>1)
 out.passed=false; out.status='visit_assignment_not_met'; out.validation_level='proposal';
end
end
