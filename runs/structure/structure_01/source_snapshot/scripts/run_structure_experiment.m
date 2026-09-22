function folder=run_structure_experiment(label,seconds,workers,windows,arms,sourceFiles)
%RUN_STRUCTURE_EXPERIMENT Paired expensive-window replacements, no seed mixing.
if nargin<1||isempty(label),label=['structure_' datestr(now,'yyyymmdd_HHMMSS')];end
if nargin<2,seconds=480;end
if nargin<3,workers=4;end
if nargin<4,windows=4;end
if nargin<5,arms='RSDP';end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
if nargin<6||isempty(sourceFiles)
 sourceFiles={fullfile(root,'runs','joint_control','joint_compare_02','seed01_B','result.mat')};
end
if ischar(sourceFiles)||isstring(sourceFiles),sourceFiles=cellstr(sourceFiles);end
validateattributes(seconds,{'numeric'},{'scalar','finite','positive'});
validateattributes(workers,{'numeric'},{'scalar','integer','>=',1,'<=',4});
validateattributes(windows,{'numeric'},{'scalar','integer','positive'});
assert(ischar(arms)&&isrow(arms)&&all(ismember(arms,'RSDP'))&&numel(unique(arms))==numel(arms)&&ismember('R',arms));
assert(~isempty(regexp(label,'^[A-Za-z0-9_-]+$','once')));
assert(license('test','Optimization_Toolbox'),'Optimization Toolbox required.');
folder=fullfile(root,'runs','structure',label);assert(~isfolder(folder),'Use a new run label.');mkdir(folder);
cfg=ctocscreen.targetedDefaults();cfg.max_wall_s=seconds;cfg.chunk_seconds=60;
cfg.joint_evaluations=12000;cfg.minimum_gap_s=30;cfg.visit_radius_km=1;
cfg.search_visit_radius_km=.99;cfg.fragment_tolerance_km=1;
cfg.zero_velocity_window=2;cfg.restart_time_sigma_s=120;cfg.time_radius=.65;
cfg.stop_file=fullfile(folder,'STOP');cfg.seed=20370921;cfg.source_hash=ctocscreen.implementationHash();
cfg.arms=arms;cfg.windows_per_source=windows;cfg.requested_workers=workers;
sources=cell(size(sourceFiles));problems=cell(size(sourceFiles));windowRows=zeros(0,7);
for b=1:numel(sourceFiles)
 d=load(sourceFiles{b});p=d.p;
 assert(size(p.states0,1)==35,'This experiment requires the 35-target competition problem.');
 if b>1,assert(isequal(p.states0,problems{1}.states0)&&p.mu_km3_s2==problems{1}.mu_km3_s2,'Sources must use the same target data and dynamics.');end
 if isfield(d,'out'),s=d.out.best.schedule;elseif isfield(d,'elite'),s=d.elite.schedule;else,error('Expected out.best or elite.');end
 ir=ctocscreen.propagateSchedule(s,p,true);nr=ctocscreen.propagateSchedule(s,p,false);assert(ir.passed&&nr.passed);
 source=struct('schedule',s,'evaluation',nr,'independent',ir);sources{b}=source;problems{b}=p;
 pl=ctocscreen.arcPlan(s,nr);u=vecnorm(s.delta_v_km_s,2,2);candidates=zeros(0,3);
 for k=1:numel(pl.counts)-2
  if all(pl.counts(k:k+2)==1)
   % Include outgoing reconnection cost in ranking where it exists.
   candidates(end+1,:)=[k sum(u(k:min(numel(u),k+3))) sum(u(k:k+2))]; %#ok<AGROW>
  end
 end
 [~,rank]=sort(candidates(:,2),'descend');selected=[];
 for j=rank'
  k=candidates(j,1);
  if isempty(selected)||all(abs(k-selected)>=3)
   selected(end+1)=k;events=cumsum(pl.counts);ids=pl.ids(events(k):events(k+2)); %#ok<AGROW>
   windowRows(end+1,:)=[b k candidates(j,2:3) ids']; %#ok<AGROW>
   if numel(selected)==windows,break;end
  end
 end
 assert(numel(selected)==windows,'Not enough disjoint three-S1 windows; reduce windows.');
end
manifest=array2table(windowRows,'VariableNames',{'source','first_arc','core_plus_outgoing_km_s','core_km_s','target1','target2','target3'});
writetable(manifest,fullfile(folder,'windows.csv'));
products=ver;save(fullfile(folder,'design.mat'),'sources','problems','sourceFiles','cfg','manifest','products','-v7.3');
snap=fullfile(folder,'source_snapshot');mkdir(snap);copyfile(fullfile(root,'src'),fullfile(snap,'src'));copyfile(fullfile(root,'scripts'),fullfile(snap,'scripts'));
copyfile(fullfile(root,'tests'),fullfile(snap,'tests'));
copyfile(fullfile(root,'docs','STRUCTURE_EXPERIMENT.md'),fullfile(snap,'DESIGN.md'));
if workers>1&&~license('test','Distrib_Computing_Toolbox'),warning('Using identical serial algorithm.');workers=1;end
cfg.actual_workers=workers;save(fullfile(folder,'design.mat'),'cfg','-append');
if workers>1
 pool=gcp('nocreate');if ~isempty(pool)&&pool.NumWorkers~=workers,error('Existing pool size differs; close it explicitly or match worker count.');end
 if isempty(pool),pool=parpool('Processes',workers);end
end
jobs=zeros(0,2);
for w=1:height(manifest),for a=circshift(1:numel(arms),[0 -(w-1)]),jobs(end+1,:)=[w a];end,end %#ok<AGROW>
save(fullfile(folder,'design.mat'),'jobs','-append');
for first=1:workers:size(jobs,1)
 if isfile(cfg.stop_file),break;end
 last=min(first+workers-1,size(jobs,1));if workers>1,futures=parallel.FevalFuture.empty;end
 for j=first:last
  w=jobs(j,1);arm=arms(jobs(j,2));b=manifest.source(w);cc=cfg;cc.seed=cfg.seed+7919*w;
  dest=fullfile(folder,sprintf('window%02d_%s',w,arm));
  if workers>1,futures(j-first+1)=parfeval(pool,@ctocscreen.runStructureTask,1,sources{b},manifest.first_arc(w),problems{b},cc,arm,dest);
  else,ctocscreen.runStructureTask(sources{b},manifest.first_arc(w),problems{b},cc,arm,dest);end
 end
 if workers>1,for j=first:last,fetchNext(futures);end,end
 report_structure_experiment(folder);
end
report_structure_experiment(folder);fprintf('Structure results: %s\n',folder);
end
