function c = bsConfig(varargin)
%BSCONFIG V1 beam-search configuration with the four expanded switches ON.
%
%   c = bsConfig()
%   c = bsConfig('tasks',8,'workers',8,'max_wall_s',480,'beam_width',8)
%
% The four switches (relative to the historical V1 defaults):
%   (1) beam_width              1  -> 8      keep W partial solutions per layer
%   (2) construct_target_count  5  -> 35     consider every remaining target
%   (3) construct_times_s       5  -> 10     denser flight-time grid
%   (4) max_revolutions         2  -> 3      more Lambert revolution branches
%
% All of these are search approximations, not competition constraints.
% Nothing in this file modifies src/+ctocscreen.

root = fileparts(fileparts(mfilename('fullpath')));   % beamsearch -> simulation
addpath(fullfile(root,'src'));

c = ctocscreen.defaultConfig('batch');

% ---- switch (1): beam width -------------------------------------------------
c.beam_width = 8;

% ---- switch (2): every remaining target is a candidate ----------------------
c.construct_target_count = 35;

% ---- switch (3): denser flight-time grid ------------------------------------
c.construct_times_s = [1800 3600 7200 10800 14400 21600 28800 43200 57600 86400];

% ---- per-leg time budgeting (same heuristic expanded_20260921_a used) --------
% Without this cap a single leg may consume the whole horizon and the replayed
% candidate can trip the "total time exceeds horizon" check by rounding.
c.construct_time_budget_factor = 1.6;

% ---- beam time reserve ------------------------------------------------------
% 0 = verbatim constructGreedy rules. Set e.g. 10 to pull every leg 10 s off the
% horizon boundary; without it a construction whose last leg lands exactly on the
% horizon can be rejected by evaluate ("Total time exceeds horizon") purely from
% floating-point rounding, which wastes a whole construction. Both the original
% constructor and the beam show this occasionally.
c.beam_time_margin_s = 0;

% ---- switch (4): more revolution branches -----------------------------------
c.branch_policy.max_revolutions = 3;

% ---- the remaining rules, set exactly as expanded_20260921_a / campaign16 ----
c.construct_refine_count = 4;                       % fminbnd polish on the best legs
c.construction_min_altitude_after_first_km = 200;   % construction altitude floor

% ---- beam fan-out control ---------------------------------------------------
% edges per layer ~= beam_width * beam_targets_per_node * beam_times_per_target
%                    * beam_branches_per_edge
% one edge ~= one Lambert enumeration + one arc check (~1-3 ms measured).
% With the defaults below one full 35-layer beam is ~17k edges, i.e. tens of
% seconds; beam_seconds caps a single construction so it stops early instead of
% eating the whole task budget.
c.beam_targets_per_node = 35;      % = construct_target_count (all remaining targets)
c.beam_times_per_target = 0;       % unused: every proposed time is used (no subsampling)
c.beam_branches_per_edge = 0;      % unused: every Lambert branch is used
c.beam_time_mode = 'grid';         % 'grid' = ctocscreen.constructGreedy's time rule
% NOTE: beam_targets_per_node / beam_times_per_target / beam_branches_per_edge and
% beam_time_mode are kept only for diagnostics; the constructor now follows
% constructGreedy's rules verbatim and the ONLY added computation is the beam
% width (it keeps cfg.beam_width partial solutions instead of one).
c.beam_endpoint_tol_km = 1.0;    % = the formal flyby limit, i.e. no stricter
                                 % filter than constructGreedy (which applies no
                                 % explicit endpoint filter at all; the Lambert
                                 % solver already aims to ~1e-4 km)
c.beam_seconds = 90;             % wall-clock cap for one beam construction
c.beam_time_throttle = 0.9;      % wall-clock safety valve only (fan-out shrink)
% 0 = verbatim constructGreedy rules. Set e.g. 10 to pull every leg 10 s off the
% horizon boundary; without it a construction whose last leg lands exactly on the
% horizon can be rejected by evaluate ("Total time exceeds horizon") purely from
% floating-point rounding, which wastes a whole construction. Both the original
% constructor and the beam show this occasionally.
c.beam_time_margin_s = 0;

% ---- task budget and parallelism -------------------------------------------
c.tasks = 8;
c.workers = 8;
c.max_wall_s = 480;              % soft search budget per task, seconds
c.max_candidates = 100000;       % outer iteration cap (beam + mutation)
c.checkpoint_interval = 10;
c.archive_size = 32;
c.export_audit_count = 32;       % independent fixed-impulse audits at the end
c.restart_probability = 0.2;     % same schedule as campaign16 (rest = mutation+SQP)
c.refine_every = 20;             % SQP refinement cadence (0 disables)
c.refine_evaluations = 2500;
c.refine_iterations = 60;
c.geometry_seed_probability = 0.9;
c.resume = false;

% ---- user overrides (name/value) -------------------------------------------
for k = 1:2:numel(varargin)
    name = varargin{k};
    value = varargin{k+1};
    if ~ischar(name) && ~isstring(name)
        error('bsConfig:name','Option names must be text.');
    end
    name = char(name);
    if strncmp(name,'branch_',7)
        c.branch_policy.(name(8:end)) = value;
    else
        c.(name) = value;
    end
end
c.search_mode = 'beam_v1';
end
