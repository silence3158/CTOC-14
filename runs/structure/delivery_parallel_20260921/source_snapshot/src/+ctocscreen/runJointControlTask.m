function out=runJointControlTask(source,p,cfg,group,folder)
%RUNJOINTCONTROLTASK One paired arm, independent incumbent, bounded budget.
if ~isfolder(folder),mkdir(folder);end
stream=RandStream('mt19937ar','Seed',cfg.seed);start=tic;
out=struct('group',group,'config',cfg,'initial',source,'best',source,'budget_best',source,'current',source, ...
 'history',[0 source.independent.total_dv_km_s 0 0], ...
 'diagnostics',{{}},'evaluations',0,'invalid_evaluations',0,'validation_rejects',0, ...
 'status','running','elapsed_s',0,'random_state',stream.State);
iteration=0;
try
 while toc(start)<cfg.max_wall_s&&~isfile(cfg.stop_file)
  iteration=iteration+1;remaining=cfg.max_wall_s-toc(start);
  cc=cfg;cc.chunk_seconds=min(cfg.chunk_seconds,remaining);
  if group=='D'&&mod(iteration,2)==0
   % Existing discrete engine supports S1/S2 only. Exploration is not a
   % claim of arbitrary-topology joint optimization.
   cc.max_wall_s=cc.chunk_seconds;cc.seed=randi(stream,2^30);
   cc.escape_enabled=true;cc.stagnation_seconds=min(15,cc.max_wall_s/3);
   cc.restart_period_seconds=min(30,cc.max_wall_s/2);
   cc.window_seconds=min(12,cc.max_wall_s);cc.checkpoint_s=30;
   trial=ctocscreen.runTargetedSearch(out.current,p,cc,fullfile(folder,sprintf('discrete_%03d',iteration)));
   out.current=trial.current;
   if trial.best.independent.total_dv_km_s<out.best.independent.total_dv_km_s,out.best=trial.best;end
   for k=1:numel(trial.diagnostics)
    dd=trial.diagnostics{k};
    if isfield(dd,'solver'),out.evaluations=out.evaluations+dd.solver.evaluations;end
   end
   diagnostic=struct('operator','discrete','elapsed_s',trial.elapsed_s, ...
    'restarts',numel(trial.restarts),'diagnostics',{trial.diagnostics});
  else
   plan=ctocscreen.arcPlan(out.current.schedule,out.current.evaluation);
   if isfield(out.current,'locked_branch_ids'),plan.locked_branch_ids=out.current.locked_branch_ids;end
   focus=randi(stream,numel(plan.counts));
   trial=ctocscreen.refineJointControl(plan,focus,p,cc,group);
   out.evaluations=out.evaluations+trial.evaluations;
   out.invalid_evaluations=out.invalid_evaluations+trial.invalid_evaluations;
   if trial.evaluation.passed&&trial.evaluation.total_dv_km_s<out.current.independent.total_dv_km_s-1e-8
    ir=ctocscreen.propagateSchedule(trial.schedule,p,true);
    nr=ctocscreen.propagateSchedule(trial.schedule,p,false);
    if ir.passed&&nr.passed
     out.current=struct('schedule',trial.schedule,'evaluation',nr,'independent',ir, ...
      'locked_branch_ids',trial.plan.locked_branch_ids);
     if ir.total_dv_km_s<out.best.independent.total_dv_km_s,out.best=out.current;end
    else,out.validation_rejects=out.validation_rejects+1;end
   end
   diagnostic=rmfield(trial,{'plan','schedule','evaluation'});diagnostic.operator='continuous';
  end
  out.diagnostics{end+1}=diagnostic;
  if toc(start)<=cfg.max_wall_s,out.budget_best=out.best;end
  out.history(end+1,:)=[toc(start) out.best.independent.total_dv_km_s out.evaluations iteration];
  checkpoint();
 end
 if isfile(cfg.stop_file),out.status='stopped';else,out.status='completed';end
catch err
 out.status='failed';out.failure_reason=err.message;out.failure_stack=err.stack;
end
checkpoint();save(fullfile(folder,'result.mat'),'out','p','cfg','-v7.3');
fprintf('JOINT %s seed %d status %s start %.9f best %.9f elapsed %.1f\n', ...
 group,cfg.seed,out.status,source.independent.total_dv_km_s,out.best.independent.total_dv_km_s,out.elapsed_s);
 function checkpoint()
  out.elapsed_s=toc(start);out.random_state=stream.State;
  save(fullfile(folder,'checkpoint.tmp.mat'),'out','p','cfg','-v7.3');
  movefile(fullfile(folder,'checkpoint.tmp.mat'),fullfile(folder,'checkpoint.mat'),'f');
 end
end
