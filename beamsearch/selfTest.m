function ok = selfTest(varargin)
%SELFTEST Fast end-to-end smoke test of the beam search (1-3 minutes).
%
%   ok = selfTest()          % tiny beam + one reference constructGreedy
%   ok = selfTest(120)       % give it more seconds
%
% It answers the three questions the first step was meant to answer:
%   1) does the beam search run on this machine at all,
%   2) how many distinct sequences and how many edge evaluations per second,
%   3) does the normal V1 pipeline (construct -> evaluate -> archive) still work.
%
% Nothing is written into src; outputs go to runs/screening/beam_selftest_*.

if nargin >= 1 && ~isempty(varargin{1})
    budget = varargin{1};
else
    budget = 90;
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
addpath(fullfile(root,'beamsearch'));

p = ctocscreen.loadProblem();
cfg = bsConfig('tasks',1,'workers',1,'max_wall_s',budget,'beam_width',2, ...
    'beam_targets_per_node',6,'beam_times_per_target',2,'beam_branches_per_edge',1, ...
    'construct_target_count',35,'construct_times_s',[3600 7200 14400 28800 43200], ...
    'construct_refine_count',0,'beam_seconds',max(30,budget*0.5), ...
    'refine_every',0,'restart_probability',1.0,'export_audit_count',2);
cfg.branch_policy.max_revolutions = 3;
cfg.run_id = "beam_selftest_" + string(datetime('now','Format','HHmmss'));

fprintf('--- SELFTEST stage 1: beam construction (W=2, %g s) ---\n', cfg.beam_seconds);
t0 = tic;
[cands, bd] = bsBeamConstruct(p, cfg, RandStream('mt19937ar','Seed',cfg.master_seed));
dt = toc(t0);
fprintf('beam status=%s, complete candidate(s)=%d, layers reached=%d, failed_depth=%d\n', ...
    bd.status, numel(cands), sum(bd.layer_nodes > 0), bd.failed_depth);
fprintf('beam edges=%d arc_checks=%d in %.1f s (%.0f edges/s)\n', ...
    bd.lambert_calls, bd.arc_checks, dt, bd.lambert_calls/max(dt,eps));

fprintf('--- SELFTEST stage 2: normal V1 pipeline (constructGreedy + evaluate) ---\n');
stream = RandStream('mt19937ar','Seed',cfg.master_seed+1);
[g, dg] = ctocscreen.constructGreedy(p, cfg, stream);
if isempty(g)
    fprintf('constructGreedy failed: %s (leg %d)\n', dg.status, dg.failure_leg);
else
    ev = ctocscreen.evaluate(g, p, cfg, true);
    fprintf('constructGreedy -> evaluate: status=%s dv=%.9f km/s visits=%d\n', ...
        ev.status, ev.total_dv_km_s, ev.unique_visit_count);
    if ~strcmp(ev.status,'two_body_verified')
        fprintf('   reason: %s (leg %d)\n', ev.failure_reason, ev.failure_leg);
    end
end

fprintf('--- SELFTEST stage 3: short search loop with the beam operator ---\n');
res = bsSearchTask(p, cfg);
fprintf('--- SELFTEST summary ---\n');
fprintf('  iterations            : %d\n', res.summary.iterations);
fprintf('  archive size          : %d\n', res.summary.n_archive);
fprintf('  best screened         : %.9f km/s\n', res.summary.best_screened_dv_km_s);
fprintf('  best independently verified: %.9f km/s\n', res.summary.best_verified_dv_km_s);
fprintf('  outputs under runs/screening/%s/\n', char(cfg.run_id));
ok = true;
end
