function batch = run_beam_v1(runId, tasks, workers, maxWallS, beamWidth, resume, varargin)
%RUN_BEAM_V1 Scaled-up V1 beam search: the four switches ON, N parallel workers.
%
%   batch = run_beam_v1()                       % 8 tasks, 8 workers, 480 s, W=8
%   batch = run_beam_v1('beam_v1_a', 8, 8, 480, 8)
%   batch = run_beam_v1('beam_v1_a', 8, 8, 480, 8, true)   % resume same run id
%
% Your own seeds (either form):
%   run_beam_v1('beam_v1_s1', 4, 8, 480, 8, false, 'master_seed', 20260922)
%       -> tasks use master_seed + 104729*(task-1)
%   run_beam_v1('beam_v1_s2', 4, 8, 480, 8, false, 'seeds', [101 202 303 404])
%       -> one explicit seed per task (cycled if fewer seeds than tasks)
%
% Any other bsConfig option can be overridden the same way, e.g.
%   ..., 'beam_width', 16, 'refine_every', 0, 'max_wall_s', 1800
%
% After the run:
%   report_beam_v1('beam_v1_a')     % figures + csv + summary.md (same style as
%                                   % scripts/report_fragment_pilot.m)
%
% In MATLAB (from anywhere):
%   cd('D:\CTOC-14\simulation')
%   addpath('beamsearch')
%   b = run_beam_v1('beam_v1_a', 8, 8, 480, 8);
%
% The four switches relative to the historical V1 search:
%   (1) beam_width             1  -> 8      keep W partial solutions per layer
%   (2) construct_target_count 5  -> 35     every remaining target is a candidate
%   (3) construct_times_s      5  -> 10     denser flight-time grid
%   (4) max_revolutions        2  -> 3      more Lambert revolution branches
%
% Nothing here modifies src/+ctocscreen; it only calls it.
%
% Stop a running batch in a controlled way: create an empty file named STOP
% inside the task directory runs/screening/<run_id>_jobNNNN/.

if nargin < 1 || isempty(runId),    runId = 'beam_v1_a'; end
if nargin < 2 || isempty(tasks),    tasks = 8; end
if nargin < 3 || isempty(workers),  workers = 8; end
if nargin < 4 || isempty(maxWallS), maxWallS = 480; end
if nargin < 5 || isempty(beamWidth),beamWidth = 8; end
if nargin < 6 || isempty(resume),   resume = false; end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
addpath(fullfile(root,'beamsearch'));

p = ctocscreen.loadProblem();
cfg = bsConfig('tasks',tasks,'workers',workers,'max_wall_s',maxWallS, ...
    'beam_width',beamWidth,'resume',resume, varargin{:});
cfg.run_id = string(runId);

fprintf('=== BEAM V1 ===\n');
fprintf('  run id        : %s\n', char(cfg.run_id));
fprintf('  tasks/workers : %d / %d\n', tasks, workers);
fprintf('  seconds/task  : %g (soft)\n', maxWallS);
fprintf('  switch 1 beam width      : %d\n', cfg.beam_width);
fprintf('  switch 2 candidate targets: %d (construct_target_count)\n', cfg.construct_target_count);
fprintf('  switch 3 time grid       : %d values\n', numel(cfg.construct_times_s));
fprintf('  switch 4 max revolutions : %d\n', cfg.branch_policy.max_revolutions);
fprintf('  beam fan-out per node    : %d targets x %d times x %d branches\n', ...
    cfg.beam_targets_per_node, cfg.beam_times_per_target, cfg.beam_branches_per_edge);
fprintf('  beam_seconds cap         : %g\n', cfg.beam_seconds);
fprintf('  total process time       : ~%g min (search only)\n', tasks*maxWallS/workers/60);

batch = bsBatch(p, cfg);

fprintf('=== DONE ===\n');
if ~isempty(batch.export_archive)
    fprintf('  best independently verified: %.12f km/s\n', batch.export_archive(1).total_dv_km_s);
else
    fprintf('  no independently verified candidate in this batch (check per-task result.mat)\n');
end
for i = 1:numel(batch.results)
    fprintf('  task%02d: iters=%d screened=%.9f verified=%.9f\n', i, ...
        batch.results{i}.summary.iterations, ...
        batch.results{i}.summary.best_screened_dv_km_s, ...
        batch.results{i}.summary.best_verified_dv_km_s);
end
end
