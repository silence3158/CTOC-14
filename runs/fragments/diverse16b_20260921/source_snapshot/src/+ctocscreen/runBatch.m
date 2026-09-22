function batch=runBatch(problem,config,resume)
%RUNBATCH Independent seeded jobs; process workers or identical serial jobs.
% Each job has its own checkpoint. Budgets apply per job, not to the batch.
if nargin<3, resume=false; end
assert(config.tasks>=1 && fix(config.tasks)==config.tasks);
useParallel=config.workers>0;
if useParallel && ~license('test','Distrib_Computing_Toolbox')
 warning('ctocscreen:batch:serial','Parallel toolbox unavailable; using serial jobs.');
 useParallel=false;
end
results=cell(1,config.tasks);
if useParallel
 pool=gcp('nocreate');
 if isempty(pool), pool=parpool('Processes',config.workers); end
 if pool.NumWorkers~=config.workers
  error('ctocscreen:batch:pool','Existing pool size differs from config.workers.');
 end
 parfor k=1:config.tasks
  results{k}=job(problem,config,k,resume);
 end
else
 for k=1:config.tasks, results{k}=job(problem,config,k,resume); end
end
archive=[];
for k=1:numel(results)
 for j=1:numel(results{k}.archive)
  entry=results{k}.archive(j);
  % Export audit repeats propagation and recomputes the actual impulse cost.
  ev=ctocscreen.evaluate(entry.candidate,problem,config,true);
  archive=ctocscreen.updateArchive(archive,entry.candidate,ev,config.archive_size);
 end
end
batch=struct('config',config,'parallel',useParallel,'results',{results},'archive',archive);
batch.export_archive=[];
for k=1:numel(results)
 for j=1:numel(results{k}.export_archive)
  e=results{k}.export_archive(j);
  batch.export_archive=ctocscreen.updateArchive(batch.export_archive,e.candidate,e.evaluation,config.archive_size);
 end
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
folder=fullfile(root,'runs','screening',char(config.run_id));
if ~isfolder(folder), mkdir(folder); end
save(fullfile(folder,'batch.mat'),'batch','-v7.3');
end
function result=job(problem,c,k,resume)
c.master_seed=mod(c.master_seed+104729*(k-1),2^32);
c.run_id=string(c.run_id)+sprintf('_job%04d',k);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
file=fullfile(root,'runs','screening',char(c.run_id),'checkpoint.mat');
state=[];
if resume && isfile(file), state=ctocscreen.loadCheckpoint(file); end
result=ctocscreen.runSearch(problem,c,state);
end
