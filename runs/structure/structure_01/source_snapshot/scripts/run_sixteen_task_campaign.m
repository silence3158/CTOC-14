function campaign=run_sixteen_task_campaign()
%RUN_SIXTEEN_TASK_CAMPAIGN Four sequential waves, four new seeds per wave.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
parent=fullfile(root,'runs','screening','campaign16_20260921');
if ~isfolder(parent), mkdir(parent); end
source=load(fullfile(root,'runs','screening','expanded_20260921_a','batch.mat'));
problem=source.batch.results{1}.problem; base=source.batch.config;
campaign=struct('baseline_dv_km_s',source.batch.export_archive(1).total_dv_km_s, ...
 'archive',source.batch.export_archive,'rounds',{{}},'completed_rounds',0, ...
 'source_hash',ctocscreen.implementationHash(),'started_utc',char(datetime('now','TimeZone','UTC')));
for wave=1:4
 c=base; c.run_id=sprintf('campaign16_20260921_wave%02d',wave);
 c.master_seed=20300000+104729*4*(wave-1);
 c.tasks=4; c.workers=4; c.max_wall_s=480; c.max_candidates=100000;
 fprintf('WAVE_START %d/4 seeds %d %d %d %d\n',wave,c.master_seed+104729*(0:3));
 path=fullfile(root,'runs','screening',c.run_id,'batch.mat');
 if isfile(path)
  loaded=load(path); b=loaded.batch;
 else
  b=ctocscreen.runBatch(problem,c,true);
 end
 for j=1:numel(b.export_archive)
  e=b.export_archive(j);
  campaign.archive=ctocscreen.updateArchive(campaign.archive,e.candidate,e.evaluation,64);
 end
 count=sum(cellfun(@(r)r.summary.iterations,b.results));
 calls=0; passed=0;
 for j=1:4
  state=ctocscreen.loadCheckpoint(b.results{j}.checkpoint_file);
  for k=1:numel(state.refinements), calls=calls+state.refinements{k}.evaluations; end
  passed=passed+numel(b.results{j}.export_archive);
 end
 summary=struct('wave',wave,'batch_path',path,'config',c,'iterations',count, ...
  'sqp_evaluations',calls,'verified_count_before_merge',passed, ...
  'best_screened_dv',b.archive(1).total_dv_km_s,'best_verified_dv',Inf);
 if ~isempty(b.export_archive), summary.best_verified_dv=b.export_archive(1).total_dv_km_s; end
 campaign.rounds{wave}=summary; campaign.completed_rounds=wave;
 elite=campaign.archive(1); config=c;
 save(fullfile(parent,'campaign.mat'),'campaign','problem','-v7.3');
 save(fullfile(parent,'elite.mat'),'elite','problem','config');
 fprintf('WAVE_DONE %d iterations%d screened%.12f verified%.12f cumulative_best%.12f\n', ...
  wave,count,summary.best_screened_dv,summary.best_verified_dv,elite.total_dv_km_s);
end
campaign.finished_utc=char(datetime('now','TimeZone','UTC'));
save(fullfile(parent,'campaign.mat'),'campaign','problem','-v7.3');
fprintf('CAMPAIGN_DONE 16_TASKS BEST %.12f\n',campaign.archive(1).total_dv_km_s);
end