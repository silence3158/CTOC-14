function batch=run_fragment_pilot(seconds,tasks,label,s2Probability)
%RUN_FRAGMENT_PILOT V2 pilot, independent tasks, no nested solver parallelism.
if nargin<1, seconds=480; end
if nargin<2, tasks=4; end
if nargin<3, label='pilot8min_20260921'; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
src=load(fullfile(root,'runs/screening/campaign16_20260921/elite.mat'));
p=src.problem; source=src.elite; cfg=ctocscreen.fragmentDefaults(); cfg.max_wall_s=seconds;
if nargin>=4, cfg.s2_probability=s2Probability; end
folder=fullfile(root,'runs/fragments',label); if ~isfolder(folder), mkdir(folder); end
cfg.source_file='runs/screening/campaign16_20260921/elite.mat';
assert(license('test','Optimization_Toolbox'),'Optimization Toolbox required.');
products=ver; save(fullfile(folder,'config.mat'),'cfg','products','p');
batch=struct('config',cfg,'baseline_dv_km_s',source.total_dv_km_s,'results',{{}},'folder',folder);
useParallel=tasks>1&&license('test','Distrib_Computing_Toolbox');
if useParallel
 pool=gcp('nocreate'); if isempty(pool), pool=parpool('Processes',min(4,tasks)); end
 f=parallel.FevalFuture.empty;
 for j=1:tasks
  c=cfg; c.seed=cfg.seed+104729*(j-1);
  f(j)=parfeval(pool,@ctocscreen.runFragmentSearch,1,source,p,c,fullfile(folder,sprintf('job%02d',j)));
 end
 for j=1:tasks
  [idx,out]=fetchNext(f); batch.results{idx}=out;
  save(fullfile(folder,'batch.mat'),'batch','-v7.3');
 end
else
 for j=1:tasks
  c=cfg; c.seed=cfg.seed+104729*(j-1);
  batch.results{j}=ctocscreen.runFragmentSearch(source,p,c,fullfile(folder,sprintf('job%02d',j)));
 end
end
batch.best_verified_dv_km_s=Inf; batch.best_job=0;
for j=1:tasks
 o=batch.results{j};
 if o.independent.passed&&o.independent.total_dv_km_s<batch.best_verified_dv_km_s
  batch.best_verified_dv_km_s=o.independent.total_dv_km_s; batch.best_job=j;
 end
end
save(fullfile(folder,'batch.mat'),'batch','-v7.3');
if batch.best_job>0
 elite=batch.results{batch.best_job}; save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
end
fprintf('PILOT_DONE tasks%d BEST %.12f job%d\n',tasks,batch.best_verified_dv_km_s,batch.best_job);
end
