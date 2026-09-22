function out=runTargetedSearch(source,p,cfg,folder)
%RUNTARGETEDSEARCH Accumulate independently verified arc fragments and repairs.
if ~isfolder(folder),mkdir(folder);end
stream=RandStream('mt19937ar','Seed',cfg.seed);start=tic;lastSave=-Inf;
if isfield(cfg,'initial_random_state'),stream.State=cfg.initial_random_state;end
out=struct('best',source,'config',cfg,'source_hash',ctocscreen.implementationHash(), ...
 'history',zeros(0,11),'discovery',{{}},'diagnostics',{{}},'restarts',{{}},'elapsed_s',0);
plan=ctocscreen.arcPlan(source.schedule,source.evaluation);candidates=[];iteration=0;generation=0;
best=source.independent.total_dv_km_s;current=source;currentCost=best;restartCount=0;lastMeaningful=0;lastRestart=0;
if isfield(cfg,'initial_best')&&cfg.initial_best.independent.passed&&cfg.initial_best.independent.total_dv_km_s<best
 out.best=cfg.initial_best;best=out.best.independent.total_dv_km_s;
end
while toc(start)<cfg.max_wall_s
 if isfile(fullfile(folder,'STOP'))||(isfield(cfg,'stop_file')&&isfile(cfg.stop_file)),break;end
 if cfg.escape_enabled&&(toc(start)-max(lastMeaningful,lastRestart)>cfg.stagnation_seconds||toc(start)-lastRestart>cfg.restart_period_seconds)
  cc=cfg;cc.kick_seconds=min(cfg.kick_seconds,max(.1,cfg.max_wall_s-toc(start)));
  [kicked,klog]=ctocscreen.kickArcSeed(current,p,cc,stream);out.restarts{end+1}=klog;lastRestart=toc(start);
  if ~isempty(kicked)
   current=kicked;currentCost=current.independent.total_dv_km_s;plan=ctocscreen.arcPlan(current.schedule,current.evaluation);
   restartCount=restartCount+1;candidates=[];
   if currentCost<best,best=currentCost;out.best=current;end
   fprintf('RESTART seed%d count%d current%.9f best%.9f changed%d\n',cfg.seed,restartCount,currentCost,best,klog.order_hamming);
  end
 end
 iteration=iteration+1;proposal=plan;M=numel(plan.counts);focus=randi(stream,M);op=1;pair=[NaN NaN];
 try
  if rand(stream)<cfg.merge_probability&&M>1
   if isempty(candidates)||mod(iteration,cfg.discovery_refresh)==0
    candidates=ctocscreen.discoverArcPairs(current.schedule,current.evaluation,p,cfg,stream);
    out.discovery{end+1}=struct('generation',generation,'candidates',candidates);
   end
   if ~isempty(candidates)
    n=min(cfg.candidate_pool,size(candidates,1));index=randi(stream,n);row=candidates(index,:);candidates(index,:)=[];
    a=row(1);b=row(2);ends=cumsum(plan.counts);m=find(ends==a,1);
    if ~isempty(m)&&m<M&&plan.counts(m)==1&&plan.counts(m+1)==1
     % Insert B next to A; the displaced targets retain their relative order.
     order=[1:a b a+1:b-1 b+1:numel(plan.ids)];proposal.ids=plan.ids(order);
     proposal.counts(m)=2;proposal.counts(m+1)=[];proposal.waits(m+1)=[];proposal.reference_v(m+1,:)=[];
     if isfield(proposal,'aim_offsets_km'),proposal.aim_offsets_km(m,:)=[];end
     focus=m;op=2;pair=[plan.ids(a) plan.ids(b)];
     lo=proposal.times(a)+60;hi=p.horizon_s;if a+2<=numel(plan.ids),hi=proposal.times(a+2)-60;end
     proposal.times(a+1)=min(hi,max(lo,row(3)));
    end
   end
  end
  cc=cfg;cc.window_seconds=min(cfg.window_seconds,max(.1,cfg.max_wall_s-toc(start)));
  trial=ctocscreen.refineArcWindow(proposal,focus,p,cc);nominal=Inf;accepted=false;independent=false;globalImproved=false;
  if trial.evaluation.passed
   nominal=trial.evaluation.total_dv_km_s;
   if nominal<currentCost-1e-7
    ir=ctocscreen.propagateSchedule(trial.schedule,p,true);independent=ir.passed;
    if ir.passed
     nr=ctocscreen.propagateSchedule(trial.schedule,p,false);
     assert(nr.passed);plan=trial.plan;plan.reference_v=nr.preburn_states(:,4:6)+trial.schedule.delta_v_km_s;
     if currentCost-ir.total_dv_km_s>=cfg.meaningful_improvement_km_s,lastMeaningful=toc(start);end
     currentCost=ir.total_dv_km_s;current=struct('schedule',trial.schedule,'evaluation',nr,'independent',ir,'arc_plan',plan);
     if currentCost<best-1e-7,best=currentCost;out.best=current;globalImproved=true;end
     generation=generation+1;accepted=true;candidates=[];
     fprintf('ACCEPT seed%d iter%d operator%d current%.12f best%.12f burns%d\n',cfg.seed,iteration,op,currentCost,best,numel(plan.counts));
    end
   end
  end
  out.history(end+1,:)=[toc(start) iteration op nominal best accepted independent numel(plan.counts) currentCost restartCount globalImproved];
  out.diagnostics{end+1}=struct('pair',pair,'focus',focus,'solver',rmfield(trial,{'schedule','evaluation','plan'}));
 catch err
  out.diagnostics{end+1}=struct('failure_reason',err.message,'stack',err.stack);
 end
 if toc(start)-lastSave>=cfg.checkpoint_s,checkpoint();lastSave=toc(start);end
end
checkpoint();save(fullfile(folder,'result.mat'),'out','p','cfg','-v7.3');
fprintf('TASK_DONE seed%d dv%.12f elapsed%.1f trials%d\n',cfg.seed,best,toc(start),iteration);
 function checkpoint()
  out.elapsed_s=toc(start);out.random_state=stream.State;out.current_plan=plan;out.current=current;out.candidates=candidates;
  save(fullfile(folder,'checkpoint.tmp.mat'),'out','p','cfg','-v7.3');
  movefile(fullfile(folder,'checkpoint.tmp.mat'),fullfile(folder,'checkpoint.mat'),'f');
  fid=fopen(fullfile(folder,'progress.txt'),'w');clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
  fprintf(fid,'seed=%d elapsed=%.1f iteration=%d best=%.12f current=%.9f best_burns=%d accepts=%d restarts=%d\n', ...
   cfg.seed,out.elapsed_s,iteration,best,currentCost,numel(out.best.schedule.maneuver_times_s),generation,restartCount);
 end
end
