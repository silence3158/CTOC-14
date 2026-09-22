function study=run_advanced_algorithms(label,seconds)
%RUN_ADVANCED_ALGORITHMS Small independent component experiments, not big batch.
if nargin<1,label='algorithms_20260921';end
if nargin<2,seconds=120;end
root=fileparts(fileparts(mfilename('fullpath')));cd(root);addpath('src');
a=load('runs/fragments/pilot8min_20260921/elite.mat');source=a.elite;p=a.p;
cfg=ctocscreen.advancedDefaults();cfg.beam_seconds=seconds;cfg.joint_seconds=seconds;cfg.multi_total_seconds=seconds;
folder=fullfile('runs/fragments',label); if ~isfolder(folder),mkdir(folder);end
cfg.source='runs/fragments/pilot8min_20260921/elite.mat';cfg.source_hash=ctocscreen.implementationHash();
study=struct('config',cfg,'baseline',source.independent.total_dv_km_s,'outputs',{{}}, ...
 'component_names',{{'free_initial_joint','global_beam','multi_fragments'}});
save(fullfile(folder,'config.mat'),'cfg','p');
functions={@ctocscreen.refineFreeInitial,@ctocscreen.expandGlobalBeam,@ctocscreen.runMultiTrials};
if license('test','Distrib_Computing_Toolbox')
 pool=gcp('nocreate');if isempty(pool),pool=parpool('Processes',3);end
 futures=parallel.FevalFuture.empty;
 for j=1:3,futures(j)=parfeval(pool,functions{j},1,source,p,cfg);end
 for j=1:3
  [idx,o]=fetchNext(futures);study.outputs{idx}=o;
  save(fullfile(folder,[study.component_names{idx} '.mat']),'o','p','cfg','-v7.3');
  save(fullfile(folder,'study.mat'),'study','-v7.3');
  fprintf('COMPONENT_DONE %s\n',study.component_names{idx});
 end
else
 for j=1:3
  o=functions{j}(source,p,cfg);study.outputs{j}=o;
  save(fullfile(folder,[study.component_names{j} '.mat']),'o','p','cfg','-v7.3');
 end
end
elite=source;
o=study.outputs{1};if o.independent.passed&&o.independent.total_dv_km_s<elite.independent.total_dv_km_s,elite=o;end
o=study.outputs{2}.best_verified;if o.independent.total_dv_km_s<elite.independent.total_dv_km_s,elite=o;end
o=study.outputs{3}.best;if o.independent.total_dv_km_s<elite.independent.total_dv_km_s,elite=o;end
study.best_verified_dv=elite.independent.total_dv_km_s;
save(fullfile(folder,'study.mat'),'study','-v7.3');save(fullfile(folder,'elite.mat'),'elite','p','cfg','-v7.3');
fprintf('ADVANCED_DONE baseline%.12f best%.12f\n',study.baseline,study.best_verified_dv);
end
