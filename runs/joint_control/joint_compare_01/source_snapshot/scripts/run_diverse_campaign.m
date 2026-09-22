function study=run_diverse_campaign(label,seconds,tasks)
%RUN_DIVERSE_CAMPAIGN Distinct whole-mission seeds; never collapse waves to one elite.
if nargin<1,label='diverse16_20260921';end
if nargin<2,seconds=480;end
if nargin<3,tasks=16;end
root=fileparts(fileparts(mfilename('fullpath')));cd(root);addpath('src');
folder=fullfile(root,'runs','fragments',label);assert(~isfolder(folder));mkdir(folder);
a=load('runs/fragments/targeted16_20260921/elite.mat');base=a.elite;elite=base;p=a.p;
cfg=ctocscreen.targetedDefaults();cfg.escape_enabled=true;cfg.max_wall_s=seconds;cfg.seed=20340921;
cfg.source_hash=ctocscreen.implementationHash();cfg.stop_file=fullfile(folder,'STOP');
snapshot=fullfile(folder,'source_snapshot');mkdir(snapshot);copyfile('src',fullfile(snapshot,'src'));copyfile('scripts',fullfile(snapshot,'scripts'));copyfile('tests',fullfile(snapshot,'tests'));
products=ver;save(fullfile(folder,'config.mat'),'cfg','p','products','base');
useParallel=license('test','Distrib_Computing_Toolbox');workers=1;
if useParallel
 workers=min(4,tasks);pool=gcp('nocreate');if ~isempty(pool)&&pool.NumWorkers~=workers,delete(pool);pool=[];end
 if isempty(pool),pool=parpool('Processes',workers);end
end
sources={};seedLogs={};orders=zeros(0,numel(base.schedule.event_target_ids));
a=load('runs/screening/campaign16_20260921/campaign.mat');ar=a.campaign.archive;
% Historical alternatives supply different order basins, not 64 near-duplicate elites.
for k=1:numel(ar)
 if numel(sources)>=min(6,tasks),break;end
 ids=ar(k).candidate.order(:)';if ismember(ids,orders,'rows'),continue;end
 s=ctocscreen.importV1Candidate(ar(k).candidate,ar(k).evaluation);
 nr=ctocscreen.propagateSchedule(s,p,false);if ~nr.passed,continue;end
 ir=ctocscreen.propagateSchedule(s,p,true);if ~ir.passed,continue;end
 sources{end+1}=struct('schedule',s,'evaluation',nr,'independent',ir);orders(end+1,:)=ids;
 seedLogs{end+1}=struct('mode','distinct_historical_order','archive_index',k,'initial_cost',ir.total_dv_km_s);
 fprintf('SEED %d historical cost%.9f changed%d\n',numel(sources),ir.total_dv_km_s,sum(ids~=base.schedule.event_target_ids'));
end
save(fullfile(folder,'seeds.mat'),'sources','seedLogs','orders','p','cfg','-v7.3');
while numel(sources)<tasks
 count=min(workers,tasks-numel(sources));pending=cell(1,count);logs=cell(1,count);
 if useParallel,futures=parallel.FevalFuture.empty;end
 for j=1:count
  seed=cfg.seed+numel(sources)*10007+j*104729;mode='kick';if mod(numel(sources)+j,3)==0,mode='greedy';end
  parent=base;if mod(j,2)==0,parent=sources{mod(j+numel(sources)-1,numel(sources))+1};end
  if useParallel,futures(j)=parfeval(pool,@ctocscreen.freshMissionSeed,2,parent,p,cfg,seed,mode);
  else,[pending{j},logs{j}]=ctocscreen.freshMissionSeed(parent,p,cfg,seed,mode);end
 end
 if useParallel
  for j=1:count,[idx,s,l]=fetchNext(futures);pending{idx}=s;logs{idx}=l;end
 end
 for j=1:count
  s=pending{j};retry=0;ids=[];if ~isempty(s),ids=s.schedule.event_target_ids(:)';end
  while isempty(s)||ismember(ids,orders,'rows')
   retry=retry+1;assert(retry<=20,'Unable to obtain distinct sequence.');
   parent=sources{mod(retry+j-1,numel(sources))+1};
   [s,logs{j}]=ctocscreen.freshMissionSeed(parent,p,cfg,cfg.seed+numel(sources)*997+retry*1009,'kick');
   if ~isempty(s),ids=s.schedule.event_target_ids(:)';end
  end
  sources{end+1}=s;seedLogs{end+1}=logs{j};orders(end+1,:)=ids;
  fprintf('SEED %d %s cost%.9f changed%d\n',numel(sources),logs{j}.mode,s.independent.total_dv_km_s,sum(ids~=base.schedule.event_target_ids'));
 end
 save(fullfile(folder,'seeds.mat'),'sources','seedLogs','orders','p','cfg','-v7.3');
end
assert(size(unique(orders,'rows'),1)==tasks);
distance=zeros(tasks);for i=1:tasks,for j=1:tasks,distance(i,j)=sum(orders(i,:)~=orders(j,:));end,end
manifest=zeros(tasks,5);
for j=1:tasks
 manifest(j,:)=[j sources{j}.independent.total_dv_km_s sum(orders(j,:)~=base.schedule.event_target_ids') ...
  sqrt(mean((sources{j}.schedule.event_times_s-base.schedule.event_times_s).^2)) min(distance(j,setdiff(1:tasks,j)))];
end
writetable(array2table(manifest,'VariableNames',{'task','initial_dv_km_s','changed_positions_from_baseline','time_rms_s','nearest_order_hamming'}),fullfile(folder,'seed_manifest.csv'));
save(fullfile(folder,'seeds.mat'),'sources','seedLogs','orders','distance','p','cfg','-v7.3');
study=struct('baseline',base.independent.total_dv_km_s,'results',{{}},'config',cfg,'tasks',tasks,'manifest',manifest);
save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
for first=1:workers:tasks
 if isfile(cfg.stop_file),break;end
 last=min(tasks,first+workers-1);fprintf('DIVERSE_WAVE %d tasks%d-%d\n',ceil(first/workers),first,last);
 if useParallel,futures=parallel.FevalFuture.empty;end
 for j=first:last
  cc=cfg;cc.seed=cfg.seed+j*7919;rad=[.35 .5 .65 .8];cc.time_radius=rad(mod(j-1,4)+1);
  if useParallel,futures(j-first+1)=parfeval(pool,@ctocscreen.runTargetedSearch,1,sources{j},p,cc,fullfile(folder,sprintf('task%02d',j)));
  else,out=ctocscreen.runTargetedSearch(sources{j},p,cc,fullfile(folder,sprintf('task%02d',j)));collect(j,out);end
 end
 if useParallel
  for j=first:last,[idx,out]=fetchNext(futures);collect(first+idx-1,out);end
 end
end
study.best_verified_dv=elite.independent.total_dv_km_s;save(fullfile(folder,'study.mat'),'study','-v7.3');
fprintf('DIVERSE_DONE baseline%.12f best%.12f\n',study.baseline,study.best_verified_dv);
 function collect(id,out)
  study.results{id}=struct('dv',out.best.independent.total_dv_km_s,'elapsed_s',out.elapsed_s, ...
   'history',out.history,'config',out.config,'burns',numel(out.best.schedule.maneuver_times_s),'restarts',{out.restarts});
  if out.best.independent.total_dv_km_s<elite.independent.total_dv_km_s,elite=out.best;elite.winning_task=id;end
  save(fullfile(folder,'study.mat'),'study','-v7.3');save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
  fprintf('COLLECT task%d start%.9f best%.9f global%.9f restarts%d\n',id,sources{id}.independent.total_dv_km_s,out.best.independent.total_dv_km_s,elite.independent.total_dv_km_s,numel(out.restarts));
 end
end
