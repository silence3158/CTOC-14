function ev = evaluate(candidate, problem, config, verify)
%EVALUATE Sequentially replay a screening candidate in the two-body model.
%   ev = ctocscreen.evaluate(candidate, problem, config, verify)
%
% Candidate references and nonzero IDs select exact branches. Zero IDs with
% no reference explicitly request greedy branch selection, not sequence DP.
% The function propagates the real arrival state forward and never resets
% the inspector state to a target state.

if nargin < 3 || isempty(config)
    config = struct();
end
% The legacy verify argument is retained; all physical checks always run.
% Independent fixed-impulse verification is ctocscreen.verifyIndependent.
config = fillDefault(config, 'flyby_limit_km', 1.0);
config = fillDefault(config, 'min_altitude_km', 200.0);
config = fillDefault(config, 'branch_policy', struct( ...
    'max_revolutions', 2, 'z_grid_points', 800, 'endpoint_tol_km', 1e-3));
config = fillDefault(config, 'arc_samples', 500);
% Configuration may tighten numerical acceptance, never relax formal limits.
config.flyby_limit_km=min(1,config.flyby_limit_km);
config.min_altitude_km=max(200,config.min_altitude_km);

tStart = tic;
N = 35;
ev = emptyEvaluation(N, string(problem.model_id));
ev.planned_count = N;

try
    q=candidate.initial_q;
    if numel(q)~=6 || any(~isfinite(q)) || q(1)<problem.re_km+590 || ...
       q(1)>problem.re_km+610 || hypot(q(2),q(3))>=0.001 || q(4)<0 || q(4)>pi || ...
       ~isscalar(candidate.wait_s) || ~isfinite(candidate.wait_s) || any(~isfinite(candidate.tof_s))
        ev.status='invalid_input'; ev.failure_reason='Initial orbit or times outside model bounds.'; return;
    end
    if numel(candidate.order) ~= N || numel(candidate.tof_s) ~= N
        ev.status = 'invalid_input';
        ev.failure_reason = 'Candidate must have 35 planned targets.';
        return;
    end
    if any(sort(candidate.order) ~= (1:N))
        ev.status = 'invalid_input';
        ev.failure_reason = 'Candidate order must be a permutation of 1..35.';
        return;
    end
    if candidate.wait_s < 0 || any(candidate.tof_s <= 0)
        ev.status = 'invalid_input';
        ev.failure_reason = 'wait_s and tof_s must be positive.';
        return;
    end
    if candidate.wait_s + sum(candidate.tof_s) > problem.horizon_s
        ev.status = 'constraint_violation';
        ev.failure_reason = 'Total time exceeds problem.horizon_s.';
        return;
    end

    mu = problem.mu_km3_s2;
    re = problem.re_km;
    state = ctocscreen.initialState(candidate.initial_q, mu, re);
    ev.initial_state = state;
    ev.propagation_calls = 0;

    arcWait = ctocscreen.checkArc(state, candidate.wait_s, problem, config);
    [state, infWait] = ctocscreen.propagateTwoBody(state, candidate.wait_s, mu);
    ev.propagation_calls = ev.propagation_calls + 1;
    if ~strcmp(infWait.status, 'ok')
        ev.status = 'solver_failure';
        ev.failure_leg = 0;
        ev.failure_reason = 'Initial coast propagation failed.';
        ev.elapsed_s = toc(tStart);
        return;
    end
    if ~strcmp(arcWait.status, 'ok')
        ev.status = 'constraint_violation';
        ev.failure_leg = 0;
        ev.failure_reason = 'Initial coast arc check failed.';
        ev.elapsed_s = toc(tStart);
        return;
    end
    if arcWait.min_altitude_km < config.min_altitude_km
        ev.status = 'constraint_violation';
        ev.failure_leg = 0;
        ev.failure_reason = 'Initial coast violates minimum altitude.';
        ev.min_altitude_km = arcWait.min_altitude_km;
        ev.elapsed_s = toc(tStart);
        return;
    end

    totalDv = 0;
    minAlt = arcWait.min_altitude_km;
    maxErr = 0;
    tRel = candidate.wait_s;
    reached = 0;
    legMinAlt = nan(1,N);

    for k = 1:N
        tRel = tRel + candidate.tof_s(k);
        [targetState, infT] = ctocscreen.targetStates( ...
            problem, candidate.order(k), tRel);
        if ~strcmp(infT(1).status, 'ok')
            ev.status = 'solver_failure';
            ev.failure_leg = k;
            ev.failure_reason = 'Target ephemeris propagation failed.';
            break;
        end
        rTarget = targetState(1:3);

        if isfield(candidate, 'branch_refs')
            if numel(candidate.branch_refs) < k || isempty(candidate.branch_refs{k})
                ev.status = 'no_branch';
                ev.failure_leg = k;
                ev.failure_reason = sprintf('No preselected safe branch reference on leg %d.', k);
                break;
            end
        end

        [branches, binfo] = ctocscreen.enumerateBranches( ...
            state(1:3), rTarget, candidate.tof_s(k), mu, config.branch_policy);
        ev.lambert_calls = ev.lambert_calls + 1;
        if isempty(branches) || ~strcmp(binfo.status, 'ok')
            ev.status = 'no_branch';
            ev.failure_leg = k;
            ev.failure_reason = sprintf('No Lambert branch on leg %d.', k);
            break;
        end

        refBranch = [];
        if isfield(candidate, 'branch_refs') && ...
                numel(candidate.branch_refs) >= k && ~isempty(candidate.branch_refs{k})
            refBranch = candidate.branch_refs{k};
        end
        if isempty(refBranch) && isfield(candidate,'branch_ids') && candidate.branch_ids(k)~=0
            idx=find([branches.branch_id]==candidate.branch_ids(k),1);
            if isempty(idx)
                ev.status='no_branch'; ev.failure_leg=k; ev.failure_reason='Requested branch ID missing.'; break;
            end
            refBranch=branches(idx);
        end
        branch = selectBranchForLeg(branches, refBranch, state(4:6));
        if isempty(branch)
            ev.status = 'no_branch';
            ev.failure_leg = k;
            ev.failure_reason = sprintf('Requested branch missing on leg %d.', k);
            break;
        end

        dv = branch.v_depart - state(4:6);
        executedVelocity = state(4:6) + dv;
        [stateNext, infNext] = ctocscreen.propagateTwoBody( ...
            [state(1:3), executedVelocity], candidate.tof_s(k), mu);
        ev.propagation_calls = ev.propagation_calls + 1;
        if ~strcmp(infNext.status, 'ok')
            ev.status = 'solver_failure';
            ev.failure_leg = k;
            ev.failure_reason = sprintf('Leg %d propagation failed.', k);
            break;
        end

        endErr = norm(stateNext(1:3) - rTarget);
        arc = ctocscreen.checkArc( ...
            [state(1:3), executedVelocity], candidate.tof_s(k), problem, config);
        legMinAlt(k) = arc.min_altitude_km;

        ev.depart_times_s(k) = tRel - candidate.tof_s(k);
        ev.arrive_times_s(k) = tRel;
        ev.departure_states(k,:) = state;
        ev.arrival_states(k,:) = stateNext;
        ev.delta_v(k,:) = dv;

        if endErr > config.flyby_limit_km
            ev.status = 'constraint_violation';
            ev.failure_leg = k;
            ev.failure_reason = sprintf('Leg %d endpoint error exceeds flyby limit.', k);
            break;
        end
        if ~strcmp(arc.status,'ok') || ~(isfinite(arc.min_altitude_km) && arc.min_altitude_km >= config.min_altitude_km)
            ev.status = 'constraint_violation';
            ev.failure_leg = k;
            ev.failure_reason = sprintf('Leg %d violates minimum altitude.', k);
            break;
        end

        totalDv = totalDv + norm(dv);
        maxErr = max(maxErr, endErr);
        minAlt = min(minAlt, arc.min_altitude_km);
        reached = reached + 1;
        state = stateNext;
    end

    ev.total_dv_km_s = totalDv;
    ev.duration_s = candidate.wait_s + sum(candidate.tof_s);
    ev.reached_planned_count = reached;
    ev.min_altitude_km = minAlt;
    ev.max_endpoint_error_km = maxErr;
    ev.leg_min_altitude_km = legMinAlt;
    ev.final_state = state;
    ev.constraint_diagnostics = struct();
    ev.constraint_diagnostics.flyby_limit_km = config.flyby_limit_km;
    ev.constraint_diagnostics.min_altitude_km = config.min_altitude_km;
    ev.constraint_diagnostics.branch_policy = config.branch_policy;
    ev.constraint_diagnostics.flyby_scan_done = false;
    ev.unique_visit_count = NaN;
    ev.flyby_scan_done = false;

    if reached == N && strcmp(ev.status, 'not_run')
        ev.status = 'two_body_verified';
        ev.validation_level = 'two_body_screen';
        ev.unique_visit_count = numel(unique(candidate.order));
        ev.flyby_scan_done = true; % All 35 actual endpoint events checked.
        ev.constraint_diagnostics.flyby_scan_done = true;
    elseif strcmp(ev.status, 'not_run')
        ev.status = 'partial';
    end
    if ~strcmp(ev.status,'two_body_verified')
        ev.partial_dv_km_s=totalDv; ev.total_dv_km_s=Inf;
    end
catch err
    ev.status = 'invalid_input';
    ev.failure_reason = string(err.message);
end

ev.elapsed_s = toc(tStart);
end

% -------------------------------------------------------------------------
function ev = emptyEvaluation(N, modelId)
ev = struct();
ev.status = 'not_run';
ev.model_id = string(modelId);
ev.validation_level = 'none';
ev.planned_count = N;
ev.reached_planned_count = 0;
ev.unique_visit_count = NaN;
ev.flyby_scan_done = false;
ev.total_dv_km_s = inf;
ev.duration_s = NaN;
ev.depart_times_s = nan(1,N);
ev.arrive_times_s = nan(1,N);
ev.departure_states = nan(N,6);
ev.arrival_states = nan(N,6);
ev.delta_v = nan(N,3);
ev.initial_state = nan(1,6);
ev.final_state = nan(1,6);
ev.min_altitude_km = inf;
ev.max_endpoint_error_km = NaN;
ev.leg_min_altitude_km = nan(1,N);
ev.constraint_diagnostics = struct();
ev.failure_leg = 0;
ev.failure_reason = "";
ev.lambert_calls = 0;
ev.propagation_calls = 0;
ev.elapsed_s = 0;
end

function config = fillDefault(config, field, value)
if ~isfield(config, field) || isempty(config.(field))
    config.(field) = value;
end
end

function branch = selectBranchForLeg(branches, refBranch, vBefore)
branch = [];
if ~isempty(refBranch)
    if isfield(refBranch,'branch_id')
        j=find([branches.branch_id]==refBranch.branch_id,1);
        if ~isempty(j) && strcmp(branches(j).geometry,refBranch.geometry) && ...
                branches(j).rev_estimate==refBranch.rev_estimate
            branch=branches(j);
        end
        return;
    end
    match = false(1, numel(branches));
    for j = 1:numel(branches)
        match(j) = strcmp(branches(j).geometry, refBranch.geometry) ...
            && branches(j).rev_estimate == refBranch.rev_estimate;
    end
    idx = find(match);
    if ~isempty(idx)
        [~, jj] = min(abs([branches(idx).z] - refBranch.z));
        branch = branches(idx(jj));
        return;
    end
    return;
end
costs = inf(1, numel(branches));
for j = 1:numel(branches)
    costs(j) = norm(branches(j).v_depart - vBefore);
end
[~, jBest] = min(costs);
branch = branches(jBest);
end


