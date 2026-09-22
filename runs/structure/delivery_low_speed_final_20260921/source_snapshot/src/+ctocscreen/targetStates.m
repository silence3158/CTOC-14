function [states, info] = targetStates(problem, ids, timesSeconds)
%TARGETSTATES Propagate target Cartesian states to requested times.
%   [states, info] = ctocscreen.targetStates(problem, ids, timesSeconds)
%   ids is 1..35 in the problem's target-index space.  timesSeconds is
%   relative to problem.epoch_utc.

ids = ids(:);
n = numel(ids);
if isscalar(timesSeconds)
    times = repmat(timesSeconds, n, 1);
else
    times = timesSeconds(:);
end
if numel(times) ~= n
    error('ctocscreen:targetStates:sizeMismatch', ...
        'ids and timesSeconds must have the same number of elements.');
end
if any(ids < 1) || any(ids > numel(problem.target_ids)) ...
        || any(ids ~= floor(ids))
    error('ctocscreen:targetStates:badIds', ...
        'Target ids must be integers in 1..N.');
end

states = nan(n, 6);
info = repmat(struct('status', 'ok', 'chi', nan, 'z', nan, ...
    'iterations', 0, 'f', nan, 'g', nan, 'fdot', nan, 'gdot', nan), n, 1);
for k = 1:n
    id = ids(k);
    [states(k,:), infk] = ctocscreen.propagateTwoBody( ...
        problem.states0(id,:), times(k), problem.mu_km3_s2);
    info(k) = infk;
    if ~strcmp(infk.status, 'ok')
        info(k).status = 'solver_failure';
    end
end
end

