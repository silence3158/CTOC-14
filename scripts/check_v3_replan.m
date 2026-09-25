function report=check_v3_replan(label)
%CHECK_V3_REPLAN One real proposal/recovery check, never a performance batch.
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
assert(~isempty(regexp(label,'^[A-Za-z0-9_-]+$','once')),'ctocscreen:v3:label','Invalid run label.');
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
folder=fullfile(sim,'runs','v3','development',label);
assert(~isfolder(folder),'ctocscreen:v3:folder','Use a new output folder.'); mkdir(folder);
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs','v3','preprocessing', ...
 'targets_20260923_release','target_ephemeris.mat'));
source=fullfile(sim,'runs','v3','search','v3_core_480_02','elite.mat'); loaded=load(source,'elite');
c=ctocscreen.v3Defaults(struct('seed',101,'replan_seconds',8,'replan_cut_count',2, ...
 'completion_steps',2,'completion_beam_width',3,'completion_outputs',2, ...
 'guided_attempts',2,'guided_branches',2,'time_grid_count',3,'time_refine_iterations',1, ...
 'joint_seconds',1,'joint_iterations',1));
stream=RandStream('mt19937ar','Seed',c.seed); pheromone=containers.Map('KeyType','char','ValueType','double');
signature=ctocscreen.v3ImplementationSignature(eph); clock=tic;
[proposals,replan]=ctocscreen.v3Replan(loaded.elite.schedule,eph,c,stream,pheromone,c.replan_seconds);
candidate=[]; diagnostic=[]; verification=[];
report=struct('source_file',source,'dataset_kind','competition_targets_nominal_model', ...
 'proposal_count',numel(proposals),'joint_called',false,'result','no_new_proposal', ...
 'elapsed_s',0,'replan_elapsed_s',replan.elapsed_s,'replan_overrun_s',replan.budget_overrun_s, ...
 'parent_dv_km_s',loaded.elite.verification.total_dv_km_s,'new_dv_km_s',NaN, ...
 'visit_count',0,'official_alignment_verified',false);
if ~isempty(proposals)
 report.joint_called=true;
 [candidate,diagnostic]=ctocscreen.v3JointOptimize(proposals{1}.schedule,eph,c);
 verification=ctocscreen.v3Verify(candidate,eph,c);
 report.result='new_proposal_not_independently_feasible';
 report.visit_count=verification.visit_count;
 if verification.passed
  report.new_dv_km_s=verification.total_dv_km_s;
  report.result='independent_nominal_j2_candidate';
  if report.new_dv_km_s<report.parent_dv_km_s-1e-9, report.result='independent_nominal_j2_improvement'; end
 end
end
report.elapsed_s=toc(clock);
save(fullfile(folder,'replan_check.mat'),'report','replan','proposals','candidate','diagnostic','verification','c','signature','-v7.3');
disp(report); fprintf('Probe %.3f s; prefix extraction %.3f s; reason: %s\n', ...
 replan.replay_seconds,replan.prefix_seconds,replan.reason);
end
