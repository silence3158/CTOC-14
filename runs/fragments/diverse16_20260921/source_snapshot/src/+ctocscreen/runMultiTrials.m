function out=runMultiTrials(source,p,cfg)
%RUNMULTITRIALS Compare generic fragment topologies on real mission windows.
stream=RandStream('mt19937ar','Seed',cfg.seed); start=tic;
out=struct('config',cfg,'trials',{{}},'fragments',{{}},'best',source, ...
 'statistics',zeros(numel(cfg.templates),5),'source_hash',ctocscreen.implementationHash());
iteration=0; baseline=source.independent.total_dv_km_s;
while toc(start)<cfg.multi_total_seconds
 iteration=iteration+1; type=mod(iteration-1,numel(cfg.templates))+1; counts=cfg.templates{type};
 N=sum(counts); k=randi(stream,numel(source.schedule.event_target_ids)-N+1);
 [xin,tin,ids,z,finish]=ctocscreen.multiSeed(source.schedule,source.evaluation,k,counts,p,stream);
 cc=cfg; cc.multi_seconds=min(cc.multi_seconds,cfg.multi_total_seconds-toc(start));
 [f,log]=ctocscreen.refineMultiFragment(xin,tin,ids,counts,z,p,cc,finish);
 st=out.statistics;st(type,1)=st(type,1)+1;st(type,5)=st(type,5)+log.evaluations;
 if f.passed
  st(type,2)=st(type,2)+1;out.fragments{end+1}=f;
  [s,r]=ctocscreen.appendFragmentSuffix(source.schedule,k-1,f,p,cfg);
  if r.passed
   st(type,3)=st(type,3)+1;
   if r.total_dv_km_s<baseline
    ir=ctocscreen.propagateSchedule(s,p,true);
    if ir.passed
     st(type,4)=st(type,4)+1; baseline=ir.total_dv_km_s;
     out.best=struct('schedule',s,'evaluation',r,'independent',ir,'fragment',f);
    end
   end
  end
 end
 out.statistics=st;out.trials{end+1}=struct('type',type,'k',k,'ids',ids,'solver',log,'passed',f.passed);
end
out.elapsed_s=toc(start);out.random_state=stream.State;out.best_verified_dv=baseline;
end
