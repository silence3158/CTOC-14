function report=check_v3_diversity(label)
%CHECK_V3_DIVERSITY Short real seed-ingestion check, no optimization iteration.
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
assert(~isempty(regexp(label,'^[A-Za-z0-9_-]+$','once')),'ctocscreen:v3:label','Invalid run label.');
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs','v3','preprocessing', ...
 'targets_20260923_release','target_ephemeris.mat'));
source=fullfile(sim,'runs','v3','search','v3_core_480_02','elite.mat');
loaded=load(source,'elite'); seed=loaded.elite.schedule;
wrapped=seed; wrapped.initial_q(6)=wrapped.initial_q(6)+2*pi;
c=ctocscreen.v3Defaults(struct('budget_s',30,'root_count',1,'colony_count',1, ...
 'initial_candidates',{{seed,wrapped}},'stop_file',[mfilename('fullpath') '.m']));
folder=fullfile(sim,'runs','v3','search',label);
state=ctocscreen.v3Search(folder,eph,c);
assert(state.iteration==0&&numel(state.initial_seed_checks)==2,'ctocscreen:v3:check','Seed checks incomplete.');
assert(all(cellfun(@(x)x.passed,state.initial_seed_checks)),'ctocscreen:v3:check','Independent check failed.');
assert(numel(state.work_pool)==1&&numel(state.reinforcement_ledger)==1,'ctocscreen:v3:check','Duplicate family.');
assert(state.feedback_history{2}.amount_per_key==0,'ctocscreen:v3:check','Repeated reward.');
report=struct('source_file',source,'folder',folder,'elapsed_s',state.elapsed_s, ...
 'budget_overrun_s',state.budget_overrun_s,'visit_count',state.elite.verification.visit_count, ...
 'total_dv_km_s',state.elite.verification.total_dv_km_s,'work_count',numel(state.work_pool), ...
 'family_count',numel(state.reinforcement_ledger),'feedback',{state.feedback_history}, ...
 'result','retained_feasible_seed','official_alignment_verified',false);
save(fullfile(folder,'diversity_check.mat'),'report'); disp(report);
end
