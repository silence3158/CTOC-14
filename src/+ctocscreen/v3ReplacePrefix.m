function beam=v3ReplacePrefix(beam,old,repaired)
%V3REPLACEPREFIX Replace an old physical prefix and its stale descendants.
drop=false(size(beam)); before=old.maneuver_times_s<old.duration_s;
for k=1:numel(beam)
 s=beam{k}.schedule; kept=s.maneuver_times_s<old.duration_s;
 drop(k)=s.duration_s>=old.duration_s&&isequal(s.initial_q,old.initial_q) ...
  &&isequal(s.maneuver_times_s(kept),old.maneuver_times_s(before)) ...
  &&isequal(s.delta_v_km_s(kept,:),old.delta_v_km_s(before,:));
end
beam=[{repaired},beam(~drop)];
end
