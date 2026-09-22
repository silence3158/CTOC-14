function study=continue_diverse_campaign(label)
%CONTINUE_DIVERSE_CAMPAIGN Explicit code-version migration; top up unfinished jobs.
root=fileparts(fileparts(mfilename('fullpath')));cd(root);addpath('src');
folder=fullfile(root,'runs','fragments',label);a=load(fullfile(folder,'study.mat'));study=a.study;
a=load(fullfile(folder,'elite.mat'));elite=a.elite;p=a.p;
seedData=load(fullfile(folder,'seeds.mat'));sources=seedData.sources;tasks=numel(sources);
cfg=ctocscreen.targetedDefaults();cfg.escape_enabled=true;cfg.max_wall_s=480;cfg.seed=20340921;
cfg.source_hash=ctocscreen.implementationHash();cfg.stop_file=fullfile(folder,'STOP_CONTINUATION');
snapshot=fullfile(folder,'source_snapshot_time_fix');assert(~isfolder(snapshot));mkdir(snapshot);
copyfile('src',fullfile(snapshot,'src'));copyfile('scripts',fullfile(snapshot,'scripts'));copyfile('tests',fullfile(snapshot,'tests'));
save(fullfile(folder,'continuation_config.mat'),'cfg','p');
pending=[];previous=cell(1,tasks);
for j=1:tasks
 path=fullfile(folder,sprintf('task%02d',j),'result.mat');
 if isfile(path)
  a=load(path);previous{j}=a.out;
  if a.out.elapsed_s>=480,continue;end
  copyfile(path,fullfile(folder,sprintf('task%02d',j),'phase1_result.mat'));
 end
 pending(end+1)=j; %#ok<AGROW>
end
useParallel=license('test','Distrib_Computing_Toolbox');workers=1;
if useParallel,workers=4;pool=gcp('nocreate');if isempty(pool),pool=parpool('Processes',workers);end,end
for offset=1:workers:numel(pending)
 ids=pending(offset:min(numel(pending),offset+workers-1));
 if isfile(cfg.stop_file),break;end
 if useParallel,futures=parallel.FevalFuture.empty;end
 for k=1:numel(ids)
  j=ids(k);cc=cfg;cc.seed=cfg.seed+j*7919;rad=[.35 .5 .65 .8];cc.time_radius=rad(mod(j-1,4)+1);
  source=sources{j};jobFolder=fullfile(folder,sprintf('task%02d',j));
  if ~isempty(previous{j})
   old=previous{j};source=old.current;cc.initial_best=old.best;cc.initial_random_state=old.random_state;
   cc.max_wall_s=max(.1,480-old.elapsed_s);jobFolder=fullfile(folder,sprintf('task%02d_continued',j));
  end
  fprintf('CONTINUE_TASK %d remaining%.3f\n',j,cc.max_wall_s);
  if useParallel,futures(k)=parfeval(pool,@ctocscreen.runTargetedSearch,1,source,p,cc,jobFolder);
  else,out=ctocscreen.runTargetedSearch(source,p,cc,jobFolder);collect(j,out);end
 end
 if useParallel
  for k=1:numel(ids),[idx,out]=fetchNext(futures);collect(ids(idx),out);end
 end
end
study.best_verified_dv=elite.independent.total_dv_km_s;study.continuation_config=cfg;
save(fullfile(folder,'study.mat'),'study','-v7.3');save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
fprintf('CONTINUATION_DONE best%.12f\n',study.best_verified_dv);
 function collect(j,out)
  if ~isempty(previous{j})
   old=previous{j};hh=out.history;hh(:,1)=hh(:,1)+old.elapsed_s;
   if ~isempty(old.history),hh(:,2)=hh(:,2)+old.history(end,2);hh(:,10)=hh(:,10)+old.history(end,10);end
   out.history=[old.history;hh];out.elapsed_s=old.elapsed_s+out.elapsed_s;
   out.diagnostics=[old.diagnostics out.diagnostics];out.restarts=[old.restarts out.restarts];out.discovery=[old.discovery out.discovery];
   out.phases={struct('source_hash',old.source_hash,'elapsed_s',old.elapsed_s,'config',old.config), ...
    struct('source_hash',out.source_hash,'elapsed_s',out.elapsed_s-old.elapsed_s,'config',out.config)};
   save(fullfile(folder,sprintf('task%02d',j),'result.mat'),'out','p','cfg','-v7.3');
  end
  study.results{j}=struct('dv',out.best.independent.total_dv_km_s,'elapsed_s',out.elapsed_s, ...
   'history',out.history,'config',out.config,'burns',numel(out.best.schedule.maneuver_times_s),'restarts',{out.restarts});
  if out.best.independent.total_dv_km_s<elite.independent.total_dv_km_s,elite=out.best;elite.winning_task=j;end
  save(fullfile(folder,'study.mat'),'study','-v7.3');save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
  fprintf('COLLECT task%d best%.9f global%.9f combined_seconds%.1f\n',j,out.best.independent.total_dv_km_s,elite.independent.total_dv_km_s,out.elapsed_s);
 end
end
