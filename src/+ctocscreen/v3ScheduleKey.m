function key=v3ScheduleKey(s)
%V3SCHEDULEKEY Exact physical-input identity; witness labels are not dynamics.
key=mat2str([s.initial_q(:);s.duration_s;s.maneuver_times_s(:);s.delta_v_km_s(:)],17);
end
