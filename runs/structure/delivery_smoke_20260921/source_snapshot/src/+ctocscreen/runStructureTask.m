function out=runStructureTask(source,first,p,cfg,arm,folder)
%RUNSTRUCTURETASK Equal-budget structural arm; fallback is not a replacement.
if ~isfolder(folder),mkdir(folder);end
switch arm
 case 'R',templates={[1 1 1]};
 case 'S',templates={[1 2],[2 1],[3]};
 case 'D',templates={[0 3]};
 case 'P',templates={[0 1 1 1]};
 otherwise,error('Unknown arm.');
end
stream=RandStream('mt19937ar','Seed',cfg.seed);start=tic;iteration=0;
out=struct('initial',source,'best',source,'budget_best',source,'best_replacement',[], ...
 'budget_replacement',[],'arm',arm,'first_arc',first,'config',cfg,'status','running', ...
 'diagnostics',{{}},'history',[0 source.independent.total_dv_km_s 0], ...
 'feasible_replacements',0,'validation_rejects',0,'elapsed_s',0);
warm=cell(size(templates));costs=inf(size(templates));
try
 while toc(start)<cfg.max_wall_s&&~isfile(cfg.stop_file)
  iteration=iteration+1;j=mod(iteration-1,numel(templates))+1;
  [pl,core]=ctocscreen.structureTemplate(source,first,templates{j},p,cfg);
  if ~isempty(warm{j}),pl=warm{j};
  elseif iteration>numel(templates)
   ix=core.first_event:core.last_event;pl.times(ix)=pl.times(ix)+cfg.restart_time_sigma_s*randn(stream,numel(ix),1);
   zero=find(pl.counts==0);pl.zero_delta_v(zero,:)=pl.zero_delta_v(zero,:)+.05*randn(stream,numel(zero),3);
   if isfield(pl,'locked_branch_ids'),pl=rmfield(pl,'locked_branch_ids');end
   [~,r]=ctocscreen.rebuildArcPlan(pl,p,cfg);if r.completed,pl.locked_branch_ids=r.branch_ids;end
  end
  cc=cfg;cc.chunk_seconds=min(cfg.chunk_seconds,cfg.max_wall_s-toc(start));
  trial=ctocscreen.refineStructureWindow(pl,core,p,cc);verified=false;
  if trial.evaluation.passed
   ir=ctocscreen.propagateSchedule(trial.schedule,p,true);nr=ctocscreen.propagateSchedule(trial.schedule,p,false);
   if ir.passed&&nr.passed
    verified=true;out.feasible_replacements=out.feasible_replacements+1;
    candidate=struct('schedule',trial.schedule,'evaluation',nr,'independent',ir, ...
     'structure_plan',trial.plan,'core',core,'template',templates{j},'cost_parts_km_s',trial.cost_parts_km_s);
    dv=ir.total_dv_km_s;
    if dv<costs(j),costs(j)=dv;warm{j}=trial.plan;end
    if isempty(out.best_replacement)||dv<out.best_replacement.independent.total_dv_km_s,out.best_replacement=candidate;end
    if dv<out.best.independent.total_dv_km_s,out.best=candidate;end
    if toc(start)<=cfg.max_wall_s
     out.budget_best=out.best;out.budget_replacement=out.best_replacement;
    end
   else,out.validation_rejects=out.validation_rejects+1;end
  end
  dd=rmfield(trial,{'plan','schedule','evaluation'});dd.template=templates{j};dd.independent_passed=verified;
  out.diagnostics{end+1}=dd;
  out.history(end+1,:)=[toc(start) out.best.independent.total_dv_km_s out.feasible_replacements];
  checkpoint();
 end
 if isfile(cfg.stop_file),out.status='stopped';else,out.status='completed';end
catch err,out.status='failed';out.failure_reason=err.message;out.failure_stack=err.stack;end
checkpoint();save(fullfile(folder,'result.mat'),'out','p','cfg','-v7.3');
fprintf('STRUCTURE %s window%d %s best %.9f replacements%d elapsed%.1f\n',arm,first,out.status,out.best.independent.total_dv_km_s,out.feasible_replacements,out.elapsed_s);
 function checkpoint()
  out.elapsed_s=toc(start);out.random_state=stream.State;out.warm_plans=warm;
  save(fullfile(folder,'checkpoint.tmp.mat'),'out','p','cfg','-v7.3');
  movefile(fullfile(folder,'checkpoint.tmp.mat'),fullfile(folder,'checkpoint.mat'),'f');
 end
end
