function [cands, diag] = bsBeamConstruct(p, cfg, stream)
%BSBEAMCONSTRUCT Width-W beam construction using exactly constructGreedy's rules.
%
%   [cands, diag] = bsBeamConstruct(problem, config, stream)
%
% Rule-for-rule the same as ctocscreen.constructGreedy:
%   * candidate targets      = the same restricted set, then every remaining
%                              target on the second pass
%   * flight-time proposals  = cfg.construct_times_s with 0.2 lognormal jitter,
%                              clamped by BOTH the remaining horizon AND
%                              cfg.construct_time_budget_factor*(T-t)/(36-k)
%   * branches               = every branch enumerateBranches returns
%   * arc acceptance         = checkArc ok AND altitude >= 200 km (after the
%                              first leg: >= cfg.construction_min_altitude_after_first_km)
%   * endpoint acceptance    = |r_end - r_target| <= cfg.beam_endpoint_tol_km
%   * continuous time polish = fminbnd on the best cfg.construct_refine_count
%                              legs, same bounds/TolX/MaxFunEvals as the original
%
% The ONLY thing that changes (and the only extra computation) is the last step:
% instead of keeping the single cheapest child, the layer keeps cfg.beam_width
% partial solutions, chosen by ctocscreen.selectDiverseBeam.
%
% Everything below adds no new numerical kernel; it only calls ctocscreen.*.
% src/+ctocscreen is never modified.

N = numel(p.target_ids);
diag = struct('status','ok','lambert_calls',0,'arc_checks',0,'children',0, ...
    'failed_depth',0,'layer_nodes',zeros(1,N),'elapsed_s',0, ...
    'refined_legs',0,'refine_improved',0);
t0 = tic;

% ---- root (same as constructGreedy) ----------------------------------------
base = ctocscreen.constructCandidates(p, struct(), stream, 1);
seedTarget = [];
if cfg.geometry_seed_probability > 0 && rand(stream) < cfg.geometry_seed_probability
    [base, seedTarget] = ctocscreen.geometrySeed(base, p, stream);
end
state = ctocscreen.initialState(base.initial_q, p.mu_km3_s2, p.re_km);
[state, info] = ctocscreen.propagateTwoBody(state, base.wait_s, p.mu_km3_s2);
if ~strcmp(info.status,'ok')
    diag.status = 'root_failure';
    cands = {};
    return;
end
root = struct('state',state,'time',base.wait_s,'order',zeros(1,0),'tof',zeros(1,0), ...
    'branch_ids',zeros(1,0),'mask',false(1,N),'cost',0);
nodes = {root};

% ---- grow the beam ---------------------------------------------------------
for depth = 1:N
    if toc(t0) >= cfg.beam_seconds
        diag.status = 'time_limit';
        break;
    end
    frac = toc(t0)/max(cfg.beam_seconds,eps);
    nTargets = cfg.construct_target_count;
    if frac > cfg.beam_time_throttle
        nTargets = max(4, round(nTargets*0.5));   % wall-clock safety valve only
    end

    children = {};
    for a = 1:numel(nodes)
        node = nodes{a};
        remaining = find(~node.mask);
        if isempty(remaining)
            continue;
        end

        % ---- time ceiling: exactly the original two rules -------------------
        maxTime = p.horizon_s - node.time - (N-depth)*cfg.min_tof_s;
        if isfield(cfg,'construct_time_budget_factor') && ~isempty(cfg.construct_time_budget_factor)
            maxTime = min(maxTime, cfg.construct_time_budget_factor ...
                *(p.horizon_s - node.time)/(36 - depth));
        end
        maxTime = maxTime - cfg.beam_time_margin_s;
        if maxTime <= cfg.min_tof_s
            continue;
        end

        floorAlt = 200;
        if depth > 1
            floorAlt = max(200, cfg.construction_min_altitude_after_first_km);
        end

        % ---- collect feasible legs (two passes, like the original) ----------
        edges = repmat(struct('id',0,'tof',0,'branch_id',0,'impulse',zeros(1,3), ...
            'state',zeros(1,6),'cost',inf), 0, 1);
        for pass = 1:2
            if pass == 1
                if numel(remaining) <= nTargets
                    tset = remaining;
                else
                    perm = randperm(stream, numel(remaining));
                    tset = remaining(perm(1:nTargets));
                    if depth == 1 && ~isempty(seedTarget) && any(remaining == seedTarget)
                        tset = unique([seedTarget, tset(1:end-1)], 'stable');
                    end
                end
            else
                tset = remaining;      % original falls back to every remaining target
            end
            for id = tset(:)'
                times = cfg.construct_times_s(:) ...
                    .* exp(0.2*randn(stream, numel(cfg.construct_times_s), 1));
                times = min(times, maxTime);
                times = unique(times(times >= cfg.min_tof_s));
                for tof = times(:)'
                    target = ctocscreen.targetStates(p, id, node.time + tof);
                    bs = ctocscreen.enumerateBranches(node.state(1:3), target(1:3), tof, ...
                        p.mu_km3_s2, cfg.branch_policy);
                    diag.lambert_calls = diag.lambert_calls + 1;
                    if isempty(bs)
                        continue;
                    end
                    for j = 1:numel(bs)
                        [edge, ok] = evalBranch(node, p, id, tof, bs(j), target, floorAlt, cfg);
                        diag.arc_checks = diag.arc_checks + 1;
                        if ok
                            edges(end+1) = edge; %#ok<AGROW>
                        end
                    end
                end
            end
            if ~isempty(edges)
                break;
            end
        end
        if isempty(edges)
            continue;
        end

        % ---- continuous time polish, same rule as constructGreedy -----------
        if cfg.construct_refine_count > 0
            [~, ord] = sort([edges.cost]);
            for t = 1:min(cfg.construct_refine_count, numel(ord))
                e = edges(ord(t));
                lo = max(cfg.min_tof_s, e.tof*0.65);
                hi = min(maxTime, e.tof*1.4);
                if hi <= lo
                    continue;
                end
                diag.refined_legs = diag.refined_legs + 1;
                [tofOpt, ~] = fminbnd(@(tt) legCost(node, p, e.id, e.branch_id, tt, floorAlt, cfg), ...
                    lo, hi, optimset('Display','off','TolX',0.5,'MaxFunEvals',24));
                [cOpt, stOpt, impOpt] = legCost(node, p, e.id, e.branch_id, tofOpt, floorAlt, cfg);
                if isfinite(cOpt) && cOpt < edges(ord(t)).cost
                    edges(ord(t)).tof = tofOpt;
                    edges(ord(t)).cost = cOpt;
                    edges(ord(t)).state = stOpt;
                    edges(ord(t)).impulse = impOpt;
                    diag.refine_improved = diag.refine_improved + 1;
                end
            end
        end

        % ---- breadth is the only added computation --------------------------
        for t = 1:numel(edges)
            e = edges(t);
            child = struct('state',e.state,'time',node.time + e.tof, ...
                'order',[node.order, e.id], 'tof',[node.tof, e.tof], ...
                'branch_ids',[node.branch_ids, e.branch_id], ...
                'mask',node.mask,'cost',node.cost + e.cost);
            child.mask(e.id) = true;
            children{end+1} = child; %#ok<AGROW>
        end
        diag.children = diag.children + numel(edges);
    end

    if isempty(children)
        diag.status = 'dead_end';
        diag.failed_depth = depth;
        break;
    end
    nodes = ctocscreen.selectDiverseBeam(children, cfg.beam_width);
    diag.layer_nodes(depth) = numel(nodes);
    fprintf('  BEAM depth%02d nodes%d children%d edges%d t%.1fs\n', ...
        depth, numel(nodes), diag.children, diag.lambert_calls, toc(t0));
end

% ---- complete candidates ---------------------------------------------------
cands = {};
for a = 1:numel(nodes)
    node = nodes{a};
    if numel(node.order) ~= N
        continue;
    end
    cands{end+1} = ctocscreen.makeCandidate(base.initial_q, node.order, base.wait_s, ...
        node.tof, node.branch_ids, base.seed, '', 'beam'); %#ok<AGROW>
end
diag.elapsed_s = toc(t0);
end

% ============================================================================
function [edge, ok] = evalBranch(node, p, id, tof, branch, target, floorAlt, cfg)
%EVALBRANCH One leg of one node: impulse, arc check, real propagation.
edge = struct('id',id,'tof',tof,'branch_id',branch.branch_id, ...
    'impulse',zeros(1,3),'state',zeros(1,6),'cost',inf);
ok = false;
impulse = branch.v_depart - node.state(4:6);
xp = [node.state(1:3), node.state(4:6) + impulse];
arc = ctocscreen.checkArc(xp, tof, p);
if ~strcmp(arc.status,'ok') || arc.min_altitude_km < floorAlt
    return;
end
[x, info] = ctocscreen.propagateTwoBody(xp, tof, p.mu_km3_s2);
if ~strcmp(info.status,'ok')
    return;
end
if norm(x(1:3) - target(1:3)) > cfg.beam_endpoint_tol_km
    return;
end
edge.impulse = impulse;
edge.state = x;
edge.cost = norm(impulse);
ok = true;
end

function [cost, stateNew, impulse] = legCost(node, p, id, branch_id, tof, floorAlt, cfg)
%LEGCOST Cost of re-flying a given branch id at a refined flight time.
cost = Inf; stateNew = zeros(1,6); impulse = zeros(1,3);
if ~(isfinite(tof) && tof >= cfg.min_tof_s)
    return;
end
target = ctocscreen.targetStates(p, id, node.time + tof);
bs = ctocscreen.enumerateBranches(node.state(1:3), target(1:3), tof, p.mu_km3_s2, cfg.branch_policy);
jj = find([bs.branch_id] == branch_id, 1);
if isempty(jj)
    return;
end
[edge, ok] = evalBranch(node, p, id, tof, bs(jj), target, floorAlt, cfg);
if ~ok
    return;
end
cost = edge.cost; stateNew = edge.state; impulse = edge.impulse;
end
