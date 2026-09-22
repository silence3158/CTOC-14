function study=run_targeted_campaign(label,seconds,tasks)
%RUN_TARGETED_CAMPAIGN Four-worker waves, sharing only independently valid elites.
if nargin<1,label='targeted16_20260921';end
if nargin<2,seconds=480;end
if nargin<3,tasks=16;end
root=fileparts(fileparts(mfilename('fullpath')));cd(root);addpath('src');
assert(license('test','Optimization_Toolbox'),'Optimization Toolbox is required.');
folder=fullfile(root,'runs','fragments',label);assert(~isfolder(folder),'Use a fresh campaign label.');mkdir(folder);
a=load('runs/fragments/algorithms_20260921/elite.mat');elite=a.elite;p=a.p;cfg=ctocscreen.targetedDefaults();cfg.max_wall_s=seconds;
cfg.source='runs/fragments/algorithms_20260921/elite.mat';cfg.source_hash=ctocscreen.implementationHash();
snapshot=fullfile(folder,'source_snapshot');mkdir(snapshot);copyfile('src',fullfile(snapshot,'src'));copyfile('scripts',fullfile(snapshot,'scripts'));copyfile('tests',fullfile(snapshot,'tests'));
products=ver;save(fullfile(folder,'config.mat'),'cfg','p','products');
study=struct('baseline',elite.independent.total_dv_km_s,'results',{{}},'config',cfg,'tasks',tasks);
parallel=license('test','Distrib_Computing_Toolbox');workers=1;
if parallel
 workers=min(4,tasks);pool=gcp('nocreate');if ~isempty(pool)&&pool.NumWorkers~=workers,delete(pool);pool=[];end
 if isempty(pool),pool=parpool('Processes',workers);end
end
for first=1:workers:tasks
 last=min(tasks,first+workers-1);source=elite;wave=ceil(first/workers);
 fprintf('WAVE_START %d tasks%d-%d baseline%.12f\n',wave,first,last,source.independent.total_dv_km_s);
 save(fullfile(folder,sprintf('wave%d_source.mat',wave)),'source','p','cfg','-v7.3');
 if parallel
  futures=parallel.FevalFuture.empty;
  for j=first:last
   cc=cfg;cc.seed=cfg.seed+j*997;cc.time_radius=[.25 .42 .65 .85];cc.time_radius=cc.time_radius(mod(j-1,4)+1);
   futures(j-first+1)=parfeval(pool,@ctocscreen.runTargetedSearch,1,source,p,cc,fullfile(folder,sprintf('task%02d',j)));
  end
  for j=first:last
   [idx,out]=fetchNext(futures);id=first+idx-1;collect(id,out);
  end
 else
  for j=first:last
   cc=cfg;cc.seed=cfg.seed+j*997;out=ctocscreen.runTargetedSearch(source,p,cc,fullfile(folder,sprintf('task%02d',j)));collect(j,out);
  end
 end
end
study.best_verified_dv=elite.independent.total_dv_km_s;save(fullfile(folder,'study.mat'),'study','-v7.3');
fprintf('CAMPAIGN_DONE baseline%.12f best%.12f\n',study.baseline,study.best_verified_dv);
 function collect(id,out)
  study.results{id}=struct('dv',out.best.independent.total_dv_km_s,'elapsed_s',out.elapsed_s, ...
   'history',out.history,'config',out.config,'burns',numel(out.best.schedule.maneuver_times_s));
  if out.best.independent.passed&&out.best.independent.total_dv_km_s<elite.independent.total_dv_km_s,elite=out.best;end
  save(fullfile(folder,'study.mat'),'study','-v7.3');save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
  fprintf('COLLECT task%d dv%.12f best%.12f\n',id,out.best.independent.total_dv_km_s,elite.independent.total_dv_km_s);
 end
end
