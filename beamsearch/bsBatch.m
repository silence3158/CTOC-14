function batch = bsBatch(p, cfg)
%BSBATCH Run cfg.tasks independent beam-search tasks, cfg.workers in parallel.
%
%   batch = bsBatch(problem, config)
%
% Every task gets its own run directory runs/screening/<run_id>_jobNNNN and its
% own seed (master_seed + 104729*(task-1), the same spacing the historical
% batches used). Results are merged at the end; only independently audited
% candidates enter the merged export archive.

root = fileparts(fileparts(mfilename('fullpath')));
tasks = cfg.tasks;
workers = cfg.workers;

useParallel = workers > 1 && license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
if workers > 1 && ~useParallel
    warning('bsBatch:noParallel', ...
        'Parallel Computing Toolbox unavailable; running the same tasks serially.');
    workers = 1;
end
if useParallel
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('Processes', workers);
    elseif pool.NumWorkers ~= workers
        error('bsBatch:pool', ...
            ['An existing pool has %d workers but %d were requested. ' ...
             'Run delete(gcp(''nocreate'')) first, then start again.'], pool.NumWorkers, workers);
    end
    fprintf('Pool ready with %d workers.\n', pool.NumWorkers);
end

% ---- per-task seeds -------------------------------------------------------
% Must be built unconditionally BEFORE any parfor: a variable that is only
% assigned inside an if-branch is not a valid broadcast variable in a parallel
% loop (MATLAB raises "Unrecognized function or variable").
taskSeeds = zeros(1,tasks);
if isfield(cfg,'seeds') && ~isempty(cfg.seeds)
    seedList = double(cfg.seeds(:)');
    for i = 1:tasks
        taskSeeds(i) = mod(seedList(mod(i-1,numel(seedList))+1), 2^32);
    end
    fprintf('Using %d explicit seed(s) supplied by the caller: %s\n', ...
        numel(seedList), mat2str(seedList));
else
    for i = 1:tasks
        taskSeeds(i) = mod(double(cfg.master_seed) + 104729*(i-1), 2^32);
    end
    fprintf('Using master seed %d with spacing 104729.\n', cfg.master_seed);
end

results = cell(1, tasks);
if useParallel
    parfor i = 1:tasks
        c = cfg;
        c.run_id = sprintf('%s_job%04d', char(cfg.run_id), i);
        c.master_seed = taskSeeds(i);
        c.tasks = 1;
        c.workers = 1;
        results{i} = bsSearchTask(p, c);
    end
else
    for i = 1:tasks
        c = cfg;
        c.run_id = sprintf('%s_job%04d', char(cfg.run_id), i);
        c.master_seed = taskSeeds(i);
        c.tasks = 1;
        c.workers = 1;
        results{i} = bsSearchTask(p, c);
    end
end

% ---- merge -----------------------------------------------------------------
parent = fullfile(root,'runs','screening',char(cfg.run_id));
if ~isfolder(parent), mkdir(parent); end
archive = repmat(struct('candidate',[],'evaluation',[],'total_dv_km_s',inf),0,1);
for i = 1:tasks
    for j = 1:numel(results{i}.export_archive)
        e = results{i}.export_archive(j);
        archive = ctocscreen.updateArchive(archive, e.candidate, e.evaluation, 64);
    end
end
batch = struct('config',cfg,'problem',p,'results',{results},'export_archive',archive, ...
    'finished_utc',char(datetime('now','TimeZone','UTC')));
save(fullfile(parent,'batch.mat'),'batch','-v7.3');
if ~isempty(archive)
    elite = archive(1);
    save(fullfile(parent,'elite.mat'),'elite','p','cfg');
end
fprintf('MERGED %d independently audited candidates into %s\n', numel(archive), parent);
end
