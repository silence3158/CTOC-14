function batch=run_expanded_search(runId,resume)
%RUN_EXPANDED_SEARCH Four broad searches with longer continuous refinement.
if nargin<1, runId='expanded_20260921'; end
if nargin<2, resume=false; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
p=ctocscreen.loadProblem(); c=ctocscreen.defaultConfig('batch');
c.run_id=string(runId); c.master_seed=20260922; c.workers=4; c.tasks=4;
c.max_candidates=500; c.max_wall_s=480; c.archive_size=32;
c.construct_target_count=35;
c.construct_times_s=[1800 3600 7200 10800 14400 21600 28800 43200 57600 86400];
c.construct_refine_count=4; c.geometry_seed_probability=.9;
c.construct_time_budget_factor=1.6;
c.construction_min_altitude_after_first_km=200;
c.branch_policy.max_revolutions=3;
c.refine_every=20; c.refine_evaluations=2500; c.refine_iterations=60;
c.restart_probability=.2; c.export_audit_count=32;
s=load(fullfile(root,'runs','screening','delivery_20260921','batch.mat'));
c.initial_candidates=arrayfun(@(e)e.candidate,s.batch.archive,'UniformOutput',false);
batch=ctocscreen.runBatch(p,c,resume);
end
