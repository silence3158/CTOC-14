function result = bsSearchTask(p, cfg)
%BSSEARCHTASK One independent beam-search task (one seed, one run directory).
%
%   result = bsSearchTask(problem, config)
%
% The loop mirrors the historical ctocscreen.runSearch structure (soft wall
% budget, STOP file, checkpoints, elite only from independent audits) but uses
% bsBeamConstruct as its global operator instead of the width-1 greedy
% constructor. src/+ctocscreen is called but never modified.

root = fileparts(fileparts(mfilename('fullpath')));
runDir = fullfile(root,'runs','screening',char(cfg.run_id));
file = fullfile(runDir,'checkpoint.mat');
if ~isfolder(runDir), mkdir(runDir); end
if isfile(file) && ~cfg.resume
    error('bsSearchTask:existingRun', ...
        'Run directory already has a checkpoint: %s (use a new run_id or resume=true).',runDir);
end

stream = RandStream('mt19937ar','Seed',cfg.master_seed);
archive = repmat(struct('candidate',[],'evaluation',[],'total_dv_km_s',inf),0,1);
history = struct('iter',0,'mode','','status','','dv',inf,'best_dv',inf, ...
    'beam_cands',0,'beam_edges',0,'elapsed_s',0);
history(1) = [];
refinements = {};
iter = 0;
elapsedBefore = 0;
signature = struct('search_mode',cfg.search_mode, ...
    'implementation_hash',ctocscreen.implementationHash(), ...
    'input_hash',p.input_hash,'model_id',p.model_id, ...
    'mu_km3_s2',p.mu_km3_s2,'re_km',p.re_km,'horizon_s',p.horizon_s);
if cfg.resume && isfile(file)
    st = ctocscreen.loadCheckpoint(file);
    archive = st.archive;
    history = st.history;
    iter = st.iter;
    stream.State = st.streamState;
    elapsedBefore = st.elapsed_s;
    if isfield(st,'refinements'), refinements = st.refinements; end
    if isfield(st,'signature') && ~isequal(st.signature.implementation_hash, signature.implementation_hash)
        warning('bsSearchTask:staleCheckpoint', ...
            'Checkpoint was written by different source code; continuing anyway.');
    end
    fprintf('[%s] resumed at iter %d, archive %d\n', char(cfg.run_id), iter, numel(archive));
end

% ---- optional warm start from previously verified elites (same as runSearch) --
if ~(cfg.resume && isfile(file)) && isfield(cfg,'initial_candidates') && ~isempty(cfg.initial_candidates)
    for j = 1:numel(cfg.initial_candidates)
        cand = cfg.initial_candidates{j};
        ev = ctocscreen.evaluate(cand, p, cfg, true);
        archive = ctocscreen.updateArchive(archive, cand, ev, cfg.archive_size);
    end
    fprintf('[%s] warm start: imported %d candidate(s), archive %d\n', ...
        char(cfg.run_id), numel(cfg.initial_candidates), numel(archive));
end

clock = tic;
while iter < cfg.max_candidates && elapsedBefore + toc(clock) < cfg.max_wall_s
    if isfile(fullfile(runDir,'STOP'))
        fprintf('[%s] STOP file found.\n', char(cfg.run_id));
        break;
    end
    iter = iter + 1;
    nBeamCands = 0;
    nBeamEdges = 0;

    if isempty(archive) || rand(stream) < cfg.restart_probability
        mode = 'beam';
        [cands, bd] = bsBeamConstruct(p, cfg, stream);
        nBeamCands = numel(cands);
        nBeamEdges = bd.lambert_calls;
        for j = 1:numel(cands)
            cands{j}.candidate_id = string(cfg.run_id) + "_beam" + iter + "_" + j;
            ev = ctocscreen.evaluate(cands{j}, p, cfg, true);
            archive = ctocscreen.updateArchive(archive, cands{j}, ev, cfg.archive_size);
        end
        status = string(bd.status);
        dv = Inf;
    else
        mode = 'mutate';
        parent = archive(randi(stream,min(4,numel(archive)))).candidate;
        cc = ctocscreen.mutateCandidate(parent, p, cfg, stream);
        [cc, dp] = ctocscreen.prepareBranches(cc, p, cfg);
        if strcmp(dp.status,'ok')
            ev = ctocscreen.evaluate(cc, p, cfg, true);
            cc.candidate_id = string(cfg.run_id) + "_mut" + iter;
            archive = ctocscreen.updateArchive(archive, cc, ev, cfg.archive_size);
            status = string(ev.status);
            dv = ev.total_dv_km_s;
        else
            status = string(dp.status);
            dv = Inf;
        end
    end

    if cfg.refine_every > 0 && mod(iter,cfg.refine_every) == 0 && ~isempty(archive)
        left = cfg.max_wall_s - elapsedBefore - toc(clock);
        if left > 5
            [rc, re, rep] = ctocscreen.refineCandidate(archive(1).candidate, p, cfg, left);
            refinements{end+1} = rep; %#ok<AGROW>
            archive = ctocscreen.updateArchive(archive, rc, re, cfg.archive_size);
            fprintf('[%s] REFINE iter%03d exitflag=%g evals=%d best=%.9f (%.1fs left)\n', ...
                char(cfg.run_id), iter, rep.exitflag, rep.evaluations, ...
                archive(1).total_dv_km_s, left);
        end
    end

    best = Inf;
    if ~isempty(archive), best = archive(1).total_dv_km_s; end
    history(end+1) = struct('iter',iter,'mode',mode,'status',status,'dv',dv, ...
        'best_dv',best,'beam_cands',nBeamCands,'beam_edges',nBeamEdges, ...
        'elapsed_s',elapsedBefore + toc(clock)); %#ok<AGROW>
    fprintf('[%s] iter%03d %-6s status=%-16s best=%.9f elapsed=%.1fs\n', ...
        char(cfg.run_id), iter, mode, char(status), best, elapsedBefore + toc(clock));
    if mod(iter,cfg.checkpoint_interval) == 0
        persist();
    end
end
persist();

% ---- independent fixed-impulse audit (export gate) -------------------------
result = struct('config',cfg,'run_dir',runDir,'checkpoint_file',file, ...
    'archive',archive,'history',history,'export_archive',[],'export_audits',{{}}, ...
    'elapsed_s',elapsedBefore + toc(clock),'summary',struct('iterations',iter, ...
    'n_archive',numel(archive),'best_screened_dv_km_s',Inf,'best_verified_dv_km_s',Inf));
if ~isempty(archive)
    result.summary.best_screened_dv_km_s = archive(1).total_dv_km_s;
end
for j = 1:min(numel(archive),cfg.export_audit_count)
    audit = ctocscreen.verifyIndependent(archive(j).candidate, p, cfg);
    result.export_audits{end+1} = audit;
    if audit.passed
        entry = archive(j);
        entry.evaluation.validation_level = 'two_body_independent_ode113';
        entry.evaluation.independent_report = audit;
        result.export_archive = [result.export_archive, entry];
    end
end
if ~isempty(result.export_archive)
    result.summary.best_verified_dv_km_s = result.export_archive(1).total_dv_km_s;
    elite = result.export_archive(1);
    save(fullfile(runDir,'elite.mat'),'elite','p','cfg');
end
save(fullfile(runDir,'result.mat'),'result');
fprintf('[%s] DONE iters=%d screened=%.9f verified=%.9f\n', char(cfg.run_id), iter, ...
    result.summary.best_screened_dv_km_s, result.summary.best_verified_dv_km_s);

    function persist()
        state = struct('iter',iter,'archive',archive,'history',history, ...
            'streamState',stream.State,'elapsed_s',elapsedBefore + toc(clock), ...
            'config',cfg,'signature',signature,'refinements',{refinements});
        ctocscreen.saveCheckpoint(file,state);
    end
end
