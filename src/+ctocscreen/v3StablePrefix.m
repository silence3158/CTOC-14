function [s,screen,report]=v3StablePrefix(s,eph,c)
%V3STABLEPREFIX Retain a fixed-pulse prefix consistent across two integrators.
% This is a numerical construction gate, not a new physical constraint.
clock=tic;
[screen,a]=ctocscreen.v3Replay(s,eph,c,false,'prefix_witnesses');
[reference,b]=ctocscreen.v3Replay(s,eph,c,true,'prefix_witnesses');
report=struct('original_duration_s',s.duration_s,'retained_duration_s',0, ...
 'max_checked_difference_km',0,'truncated',false,'passed',false, ...
 'replay_count',2,'target_reintegrations',0,'elapsed_s',0,'scope','prefix_witnesses');
safeEnd=0; safeArcs=0;
for k=1:min(numel(a.arcs),numel(b.arcs))
 sa=a.arcs{k}; sb=b.arcs{k};
 if sa.x(1)~=sb.x(1)||sa.x(end)~=sb.x(end), break; end
 tt=unique([linspace(sa.x(1),sa.x(end),max(3,ceil(diff(sa.x([1 end]))/c.scan_step_s)+1)), ...
  s.witness_times_s(s.witness_times_s>=sa.x(1)&s.witness_times_s<=sa.x(end)).']);
 xa=deval(sa,tt); xb=deval(sb,tt); dr=max(vecnorm(xa(1:3,:)-xb(1:3,:),2,1));
 report.max_checked_difference_km=max(report.max_checked_difference_km,dr);
 ha=a.heights{k}; hb=b.heights{k};
 if dr>c.prefix_consistency_km||~ha.passed||~hb.passed, break; end
 safeEnd=sa.x(end); safeArcs=k;
end
if safeEnd<=0, report.elapsed_s=toc(clock); return; end
if safeEnd<s.duration_s
 report.truncated=true; s.duration_s=safeEnd;
 keep=s.maneuver_times_s<safeEnd; s.maneuver_times_s=s.maneuver_times_s(keep);
 s.delta_v_km_s=s.delta_v_km_s(keep,:); s.witness_times_s(s.witness_times_s>safeEnd)=NaN;
 screen=retained(screen,a,safeArcs,s); reference=retained(reference,b,safeArcs,s);
end
report.passed=screen.height_passed&&reference.height_passed;
report.retained_duration_s=safeEnd;
if report.passed
 valid=screen.distance_km<=1 & reference.distance_km<=1;
 s.witness_times_s=screen.witness_times_s; s.witness_times_s(~valid)=NaN;
 screen.distance_km(~valid)=Inf; screen.visit_count=sum(valid);
 screen.passed=screen.height_passed&&all(valid);
 % Continue from the independently replayed actual state, never a target or
 % an auxiliary shooting node. Both integrations still have to agree above.
 screen.final_state=reference.final_state;
 if ~screen.passed, screen.status='constraints_not_met'; end
end
report.elapsed_s=toc(clock);
report.schedule_key=ctocscreen.v3ScheduleKey(s);
report.witness_key=mat2str(s.witness_times_s,17);
report.inclination_penalty=ctocscreen.v3TracePlanePenalty(s,b,c);
if report.passed
 [report.experience_keys,report.policy_actions]=ctocscreen.v3TraceExperience(s,b,c);
 report.experience_config=struct('pheromone_state_scales',c.pheromone_state_scales);
end
end

function out=retained(out,trace,count,s)
% Reuse the already integrated safe arcs; exclude any burn at the removed suffix.
out.distance_km=inf(35,1); out.witness_times_s=nan(35,1);
out.min_altitude_lower_km=Inf; out.sampled_min_altitude_km=Inf; out.height_passed=true;
for k=1:count
 z=trace.encounters{k}; better=z.distance_km<out.distance_km;
 out.distance_km(better)=z.distance_km(better); out.witness_times_s(better)=z.time_s(better);
 h=trace.heights{k}; out.height_passed=out.height_passed&&h.passed;
 out.min_altitude_lower_km=min(out.min_altitude_lower_km,h.lower_bound_altitude_km);
 out.sampled_min_altitude_km=min(out.sampled_min_altitude_km,h.sampled_min_altitude_km);
end
y=deval(trace.arcs{count},s.duration_s); out.final_state=y(1:6).';
out.total_dv_km_s=sum(vecnorm(s.delta_v_km_s,2,2)); out.visit_count=sum(out.distance_km<=1);
out.passed=out.height_passed&&out.visit_count==35; out.status='constraints_not_met';
if out.passed, out.status='complete'; end
out.failure_reason='';
if isfield(out,'failure_id'), out=rmfield(out,'failure_id'); end
end
