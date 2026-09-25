function value=v3RecoveryMerit(replay,ids,radius)
%V3RECOVERYMERIT Physical feasibility progress, never a fuel objective.
value=Inf;
if isempty(replay)||~replay.height_passed, return; end
distance=replay.distance_km;
if isfield(replay,'planning_distance_km'), distance=replay.planning_distance_km; end
if any(~isfinite(distance(ids))), return; end
value=norm(max(0,distance(ids)-radius));
end
