function folder=run_joint_control_experiment(label,seconds,workers,seedCount,groups)
%RUN_JOINT_CONTROL_EXPERIMENT Paired, independent A/B/C/D experiment.
% No arguments: 4 seeds x 4 arms x 480 s, max 4 process workers.
if nargin<1||isempty(label),label=['joint_' datestr(now,'yyyymmdd_HHMMSS')];end
if nargin<2,seconds=480;end
if nargin<3,workers=4;end
if nargin<4,seedCount=4;end
if nargin<5,groups='ABCD';end
validateattributes(seconds,{'numeric'},{'scalar','finite','positive'});
validateattributes(workers,{'numeric'},{'scalar','integer','positive','<=',4});
validateattributes(seedCount,{'numeric'},{'scalar','integer','positive','<=',7});
assert(ischar(groups)&&isrow(groups)&&all(ismember(groups,'ABCD'))&&numel(unique(groups))==numel(groups));
assert(~isempty(regexp(label,'^[A-Za-z0-9_-]+$','once')),'Use a simple unique label.');
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'src'));
assert(license('test','Optimization_Toolbox')&&exist('fmincon','file')==2,'Optimization Toolbox required.');
folder=fullfile(root,'runs','joint_control',label);assert(~isfolder(folder),'Use a new output label.');
mkdir(folder);
d=load(fullfile(root,'runs','fragments','diverse16b_20260921','elite.mat'));
h=load(fullfile(root,'runs','fragments','diverse16b_20260921','seeds.mat'),'sources');
p=d.p;seeds=[{d.elite} h.sources(1:seedCount-1)];
cfg=ctocscreen.targetedDefaults();cfg.max_wall_s=seconds;cfg.chunk_seconds=60;
cfg.joint_evaluations=12000;cfg.joint_iterations=150;cfg.minimum_gap_s=30;
cfg.visit_radius_km=1;cfg.search_visit_radius_km=.99;cfg.fragment_tolerance_km=1;
cfg.stop_file=fullfile(folder,'STOP');cfg.source_hash=ctocscreen.implementationHash();
cfg.scope='paired A/B/C continuous ablation; D adds existing S1/S2 discrete operators';
cfg.requested_workers=workers;cfg.groups=groups;cfg.seed_count=seedCount;
% Validate once outside the common per-arm timer. Every arm gets exactly the
% same locked-branch, independently verified starting trajectory per seed.
manifest=zeros(seedCount,10);
for k=1:seedCount
 s=seeds{k}.schedule;nr=ctocscreen.propagateSchedule(s,p,false);assert(nr.passed);
 plan=ctocscreen.arcPlan(s,nr);[ss,rr]=ctocscreen.rebuildArcPlan(plan,p,cfg);assert(rr.passed);
 assert(abs(rr.total_dv_km_s-sum(vecnorm(s.delta_v_km_s,2,2)))<1e-5,'Seed branch reconstruction changed cost.');
 ir=ctocscreen.propagateSchedule(ss,p,true);assert(ir.passed,'Seed independent verification failed.');
 nr=ctocscreen.propagateSchedule(ss,p,false);assert(nr.passed);
 seeds{k}=struct('schedule',ss,'evaluation',nr,'independent',ir,'locked_branch_ids',rr.branch_ids);
 manifest(k,:)=[k ir.total_dv_km_s numel(ss.maneuver_times_s) ss.duration_s/86400 ss.initial_q];
end
seed_manifest=array2table(manifest,'VariableNames',{'seed','initial_dv_km_s','burns','days','a_km','ec','es','i_rad','Omega_rad','u_rad'});
writetable(seed_manifest,fullfile(folder,'seed_manifest.csv'));
products=ver;save(fullfile(folder,'design.mat'),'seeds','p','cfg','products','seed_manifest','-v7.3');
snap=fullfile(folder,'source_snapshot');mkdir(snap);
copyfile(fullfile(root,'src'),fullfile(snap,'src'));copyfile(fullfile(root,'scripts'),fullfile(snap,'scripts'));
copyfile(fullfile(root,'docs','JOINT_CONTROL_EXPERIMENT.md'),fullfile(snap,'DESIGN.md'));
% Rotate group order across seed blocks. Each task starts independently;
% no arm or wave can inherit another arm's improved incumbent.
jobs=zeros(seedCount*numel(groups),2);at=0;
for k=1:seedCount
 order=circshift(1:numel(groups),[0 -(k-1)]);
 for j=order,at=at+1;jobs(at,:)=[k j];end
end
if workers>1&&~license('test','Distrib_Computing_Toolbox')
 warning('ctocscreen:serial','Parallel toolbox unavailable: same experiment runs serially.');workers=1;
end
cfg.actual_workers=workers;save(fullfile(folder,'design.mat'),'cfg','jobs','-append');
if workers>1
 pool=gcp('nocreate');
 if ~isempty(pool)&&pool.NumWorkers~=workers,error('Existing pool size differs. Close it with delete(gcp(''nocreate'')) or request its size.');end
 if isempty(pool),pool=parpool('Processes',workers);end
end
for first=1:workers:size(jobs,1)
 if isfile(cfg.stop_file),break;end
 last=min(size(jobs,1),first+workers-1);
 if workers>1,futures=parallel.FevalFuture.empty;end
 for j=first:last
  k=jobs(j,1);group=groups(jobs(j,2));cc=cfg;cc.seed=20360921+k*7919;
  taskFolder=fullfile(folder,sprintf('seed%02d_%s',k,group));
  if workers>1
   futures(j-first+1)=parfeval(pool,@ctocscreen.runJointControlTask,1,seeds{k},p,cc,group,taskFolder);
  else,ctocscreen.runJointControlTask(seeds{k},p,cc,group,taskFolder);end
 end
 if workers>1,for j=first:last,fetchNext(futures);end,end
 report_joint_control_experiment(folder);
end
report_joint_control_experiment(folder);
fprintf('Results: %s\n',folder);
end
