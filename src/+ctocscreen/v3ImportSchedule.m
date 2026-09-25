function s=v3ImportSchedule(old,m)
%V3IMPORTSCHEDULE Explicit V2 schedule -> V3 proposal; no inherited feasibility.
assert(isfield(old,'maneuver_times_s')&&isfield(old,'delta_v_km_s'), ...
 'ctocscreen:v3:import','Pass an executed V2 schedule, not a nominal visit plan.');
s=struct('schema_version','free_maneuver_v3','dynamics_id','central_j2', ...
 'initial_q',old.initial_q,'maneuver_times_s',old.maneuver_times_s, ...
 'delta_v_km_s',old.delta_v_km_s,'duration_s',old.duration_s, ...
 'witness_times_s',nan(35,1),'validation_level','proposal');
if isfield(old,'witness_times_s')
 s.witness_times_s=old.witness_times_s;
elseif isfield(old,'event_target_ids')
 for k=1:numel(old.event_target_ids)
  s.witness_times_s(old.event_target_ids(k))=old.event_times_s(k);
 end
end
s=ctocscreen.v3Normalize(s,m);
end
