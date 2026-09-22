function result=runSearch(problem,config,resumeState)
%RUNSEARCH Bounded reproducible multi-start search with local improvements.
if nargin<2, config=struct(); end
base=ctocscreen.defaultConfig(); names=fieldnames(base);
for i=1:numel(names)
 if ~isfield(config,names{i}), config.(names{i})=base.(names{i}); end
end
if nargin<3, resumeState=[]; end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
runDir=fullfile(root,'runs','screening',char(config.run_id));
file=fullfile(runDir,'checkpoint.mat');
if isempty(resumeState) && isfile(file)
 error('ctocscreen:search:existingRun','Run exists; explicitly load checkpoint to resume.');
end
if ~isfolder(runDir), mkdir(runDir); end
stream=RandStream('mt19937ar','Seed',config.master_seed);
archive=[]; history=struct('iter',{},'status',{},'failure_leg',{},'dv',{},'best_dv',{});
iter=0; elapsedBefore=0; refinements={};
if isempty(resumeState) && isfield(config,'initial_candidates')
 for j=1:numel(config.initial_candidates)
  cand=config.initial_candidates{j}; ev=ctocscreen.evaluate(cand,problem,config,true);
  archive=ctocscreen.updateArchive(archive,cand,ev,config.archive_size);
 end
end
signature=struct('input_hash',problem.input_hash,'model',problem.model_id, ...
 'implementation_hash',ctocscreen.implementationHash(), ...
 'mu',problem.mu_km3_s2,'re',problem.re_km,'states',problem.states0, ...
 'horizon',problem.horizon_s,'config',config);
if ~isempty(resumeState)
 old=resumeState.signature; current=signature;
 fields={'max_candidates','max_wall_s','display_interval','checkpoint_interval'};
 old.config=rmfield(old.config,fields); current.config=rmfield(current.config,fields);
 assert(isequaln(old,current),'ctocscreen:resume:mismatch','Checkpoint problem/config mismatch.');
 archive=resumeState.archive; history=resumeState.history;
 iter=resumeState.iter; stream.State=resumeState.streamState;
 elapsedBefore=resumeState.elapsed_s; refinements=resumeState.refinements;
end
clock=tic;
while iter<config.max_candidates && elapsedBefore+toc(clock)<config.max_wall_s
 if isfile(fullfile(runDir,'STOP')), break; end
 iter=iter+1; dp=[];
 if isempty(archive) || rand(stream)<config.restart_probability
  [candidate,diag]=ctocscreen.constructGreedy(problem,config,stream);
  if isempty(candidate)
   ev=struct('status',diag.status,'failure_leg',diag.failure_leg,'total_dv_km_s',Inf);
  else
   ev=ctocscreen.evaluate(candidate,problem,config,true);
  end
 else
  parent=archive(randi(stream,min(4,numel(archive)))).candidate;
  candidate=ctocscreen.mutateCandidate(parent,problem,config,stream);
  [candidate,dp]=ctocscreen.prepareBranches(candidate,problem,config);
  if strcmp(dp.status,'ok')
   ev=ctocscreen.evaluate(candidate,problem,config,true);
  else
   ev=struct('status',dp.status,'failure_leg',dp.failure_leg,'total_dv_km_s',Inf);
  end
 end
 if strcmp(ev.status,'two_body_verified')
  candidate.candidate_id=string(config.run_id)+"_"+iter;
  archive=ctocscreen.updateArchive(archive,candidate,ev,config.archive_size);
 end
 if config.refine_every>0 && mod(iter,config.refine_every)==0 && ~isempty(archive)
  remaining=config.max_wall_s-elapsedBefore-toc(clock);
  if remaining>0
   [rc,re,report]=ctocscreen.refineCandidate(archive(1).candidate,problem,config,remaining);
   archive=ctocscreen.updateArchive(archive,rc,re,config.archive_size);
   refinements{end+1}=report;
  end
 end
 best=Inf; if ~isempty(archive), best=archive(1).total_dv_km_s; end
 history(end+1)=struct('iter',iter,'status',ev.status,'failure_leg',ev.failure_leg, ...
  'dv',ev.total_dv_km_s,'best_dv',best);
 if mod(iter,config.display_interval)==0
  fprintf('iter %d status=%s leg=%d best=%.9f km/s elapsed=%.1fs\n', ...
   iter,ev.status,ev.failure_leg,best,elapsedBefore+toc(clock));
 end
 if mod(iter,config.checkpoint_interval)==0, persist(); end
end
persist();
result=struct('config',config,'problem',problem,'archive',archive,'history',history, ...
 'run_dir',runDir,'checkpoint_file',file,'elapsed_s',elapsedBefore+toc(clock), ...
 'summary',struct('n_archive',numel(archive),'best_dv_km_s',Inf,'iterations',iter));
if ~isempty(archive)
 result.summary.best_dv_km_s=archive(1).total_dv_km_s;
 screened_elite=archive(1); save(fullfile(runDir,'screened_elite.mat'),'screened_elite','problem','config');
end
result.export_archive=[]; result.export_audits={};
for j=1:min(numel(archive),config.export_audit_count)
 audit=ctocscreen.verifyIndependent(archive(j).candidate,problem,config);
 result.export_audits{end+1}=audit;
 if audit.passed
  entry=archive(j); entry.evaluation.validation_level='two_body_independent_ode113';
  entry.evaluation.independent_report=audit;
  result.export_archive=[result.export_archive entry];
 end
end
if ~isempty(result.export_archive)
 elite=result.export_archive(1); save(fullfile(runDir,'elite.mat'),'elite','problem','config');
elseif isfile(fullfile(runDir,'elite.mat'))
 % Do not leave a previous export presented as the latest result.
 movefile(fullfile(runDir,'elite.mat'),fullfile(runDir,'previous_elite.mat'),'f');
end
save(fullfile(runDir,'result.mat'),'result');
 function persist()
  state=struct('iter',iter,'archive',archive,'history',history,'streamState',stream.State, ...
   'config',config,'signature',signature,'elapsed_s',elapsedBefore+toc(clock),'refinements',{refinements});
  ctocscreen.saveCheckpoint(file,state);
 end
end
