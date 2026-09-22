function out=runTargetedSearch(source,p,cfg,folder)
%RUNTARGETEDSEARCH Accumulate independently verified arc fragments and repairs.
if ~isfolder(folder),mkdir(folder);end
stream=RandStream('mt19937ar','Seed',cfg.seed);start=tic;lastSave=-Inf;
out=struct('best',source,'config',cfg,'source_hash',ctocscreen.implementationHash(), ...
 'history',zeros(0,8),'discovery',{{}},'diagnostics',{{}},'elapsed_s',0);
plan=ctocscreen.arcPlan(source.schedule,source.evaluation);candidates=[];iteration=0;generation=0;
best=source.independent.total_dv_km_s;
while toc(start)<cfg.max_wall_s
 iteration=iteration+1;proposal=plan;M=numel(plan.counts);focus=randi(stream,M);op=1;pair=[NaN NaN];
 try
  if rand(stream)<cfg.merge_probability&&M>1
   if isempty(candidates)||mod(iteration,cfg.discovery_refresh)==0
    candidates=ctocscreen.discoverArcPairs(out.best.schedule,out.best.evaluation,p,cfg,stream);
    out.discovery{end+1}=struct('generation',generation,'candidates',candidates);
   end
   if ~isempty(candidates)
    n=min(cfg.candidate_pool,size(candidates,1));index=randi(stream,n);row=candidates(index,:);candidates(index,:)=[];
    a=row(1);b=row(2);ends=cumsum(plan.counts);m=find(ends==a,1);
    if ~isempty(m)&&m<M&&plan.counts(m)==1&&plan.counts(m+1)==1
     % Insert B next to A; the displaced targets retain their relative order.
     order=[1:a b a+1:b-1 b+1:numel(plan.ids)];proposal.ids=plan.ids(order);
     proposal.counts(m)=2;proposal.counts(m+1)=[];proposal.waits(m+1)=[];proposal.reference_v(m+1,:)=[];
     focus=m;op=2;pair=[plan.ids(a) plan.ids(b)];
     lo=proposal.times(a)+60;hi=p.horizon_s;if a+2<=numel(plan.ids),hi=proposal.times(a+2)-60;end
     proposal.times(a+1)=min(hi,max(lo,row(3)));
    end
   end
  end
  cc=cfg;cc.window_seconds=min(cfg.window_seconds,max(.1,cfg.max_wall_s-toc(start)));
  trial=ctocscreen.refineArcWindow(proposal,focus,p,cc);nominal=Inf;accepted=false;independent=false;
  if trial.evaluation.passed
   nominal=trial.evaluation.total_dv_km_s;
   if nominal<best-1e-7
    ir=ctocscreen.propagateSchedule(trial.schedule,p,true);independent=ir.passed;
    if ir.passed
     nr=ctocscreen.propagateSchedule(trial.schedule,p,false);
     assert(nr.passed);plan=trial.plan;plan.reference_v=nr.preburn_states(:,4:6)+trial.schedule.delta_v_km_s;
     best=ir.total_dv_km_s;out.best=struct('schedule',trial.schedule,'evaluation',nr,'independent',ir,'arc_plan',plan);
     generation=generation+1;accepted=true;candidates=[];
     fprintf('ACCEPT seed%d iter%d operator%d dv%.12f burns%d\n',cfg.seed,iteration,op,best,numel(plan.counts));
    end
   end
  end
  out.history(end+1,:)=[toc(start) iteration op nominal best accepted independent numel(plan.counts)];
  out.diagnostics{end+1}=struct('pair',pair,'focus',focus,'solver',rmfield(trial,{'schedule','evaluation','plan'}));
 catch err
  out.diagnostics{end+1}=struct('failure_reason',err.message,'stack',err.stack);
 end
 if toc(start)-lastSave>=cfg.checkpoint_s,checkpoint();lastSave=toc(start);end
end
checkpoint();save(fullfile(folder,'result.mat'),'out','p','cfg','-v7.3');
fprintf('TASK_DONE seed%d dv%.12f elapsed%.1f trials%d\n',cfg.seed,best,toc(start),iteration);
 function checkpoint()
  out.elapsed_s=toc(start);out.random_state=stream.State;out.current_plan=plan;out.candidates=candidates;
  save(fullfile(folder,'checkpoint.tmp.mat'),'out','p','cfg','-v7.3');
  movefile(fullfile(folder,'checkpoint.tmp.mat'),fullfile(folder,'checkpoint.mat'),'f');
  fid=fopen(fullfile(folder,'progress.txt'),'w');clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
  fprintf(fid,'seed=%d elapsed=%.1f iteration=%d dv=%.12f burns=%d accepts=%d\n',cfg.seed,out.elapsed_s,iteration,best,numel(plan.counts),generation);
 end
end
