function [beam,removed]=v3PruneUnstableBeam(beam,safe)
%V3PRUNEUNSTABLEBEAM Restart stored extensions of a rolled-back prefix.
% A finite beam heuristic only: different initial orbits/prefixes survive.
drop=false(size(beam)); s=safe.schedule;
for k=1:numel(beam)
 other=beam{k}.schedule;
 if other.duration_s<=s.duration_s||~isequal(other.initial_q,s.initial_q), continue; end
 before=other.maneuver_times_s<s.duration_s;
 drop(k)=isequal(other.maneuver_times_s(before),s.maneuver_times_s) ...
  && isequal(other.delta_v_km_s(before,:),s.delta_v_km_s);
end
removed=sum(drop); beam=beam(~drop);
end
