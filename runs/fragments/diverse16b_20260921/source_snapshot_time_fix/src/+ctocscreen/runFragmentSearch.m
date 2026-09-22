function out=runFragmentSearch(source,p,cfg,folder)
%RUNFRAGMENTSEARCH Reproducible pilot, S2 shots and S1 free-time local repair.
if ~isfolder(folder), mkdir(folder); end
stream=RandStream('mt19937ar','Seed',cfg.seed);
c=source.candidate; e=source.evaluation;
plan=struct('initial_q',c.initial_q,'ids',c.order(:),'times',e.arrive_times_s(:), ...
 'waits',[c.wait_s;zeros(34,1)],'skip',false(35,1), ...
 'fixed_dv',nan(35,3),'reference_v',e.departure_states(:,4:6)+e.delta_v);
s=ctocscreen.importV1Candidate(c,e); r=ctocscreen.propagateSchedule(s,p,false); assert(r.passed);
out=struct('config',cfg,'source_hash',ctocscreen.implementationHash(),'input_hash',p.input_hash, ...
 'schedule',s,'evaluation',r,'plan',plan,'history',[0 r.total_dv_km_s], ...
 'statistics',struct('iterations',0,'shots',0,'gated',0,'refined',0,'fragments',0, ...
 'reconnected',0,'s2_improvements',0,'s1_improvements',0,'solver_evaluations',0), ...
 'fragment_archive',{{}},'diagnostics',{{}},'scope',cfg.scope);
checkpoint=fullfile(folder,'checkpoint.mat'); elapsed=0;
if isfile(checkpoint)
 old=load(checkpoint,'out');
 assert(strcmp(old.out.source_hash,out.source_hash)&&strcmp(old.out.input_hash,p.input_hash));
 assert(isequal(old.out.config,cfg),'Incompatible resume configuration.');
 out=old.out; elapsed=out.search_elapsed_s; stream.State=out.random_state;
 plan=out.plan; s=out.schedule; r=out.evaluation;
end
start=tic; cpuStart=cputime; lastCheckpoint=0;
while elapsed+toc(start)<cfg.max_wall_s
 remaining=cfg.max_wall_s-elapsed-toc(start); cc=cfg;
 cc.fragment_seconds=min(cc.fragment_seconds,remaining); cc.local_seconds=min(cc.local_seconds,remaining);
 st=out.statistics; st.iterations=st.iterations+1;
 if rand(stream)<cfg.s2_probability
  [f,trial,d]=ctocscreen.shootFragments(plan,r,p,cc,stream);
  st.shots=st.shots+1; st.gated=st.gated+d.gated; st.refined=st.refined+d.refined;
  if isfield(f,'evaluations'), st.solver_evaluations=st.solver_evaluations+f.evaluations; end
  if f.passed
   st.fragments=st.fragments+1; out.fragment_archive{end+1}=struct('fragment',f,'diagnostic',d);
   [ss,rr]=ctocscreen.rebuildSchedule(trial,p,cfg);
   if rr.passed
    st.reconnected=st.reconnected+1;
    if rr.total_dv_km_s<r.total_dv_km_s
     plan=trial; s=ss; r=rr; st.s2_improvements=st.s2_improvements+1;
    end
   end
  end
  d.fragment=f; out.diagnostics{end+1}=d;
 else
  previous=r.total_dv_km_s;
  [plan,s,r,d]=ctocscreen.repairFragmentWindow(plan,s,r,p,cc,stream);
  st.solver_evaluations=st.solver_evaluations+d.evaluations;
  st.s1_improvements=st.s1_improvements+(r.total_dv_km_s<previous);
  out.diagnostics{end+1}=struct('mode','S1','solver',d);
 end
 out.statistics=st; out.plan=plan; out.schedule=s; out.evaluation=r;
 out.history(end+1,:)=[elapsed+toc(start),r.total_dv_km_s];
 if toc(start)-lastCheckpoint>=cfg.checkpoint_interval_s
  persist(); lastCheckpoint=toc(start);
  fprintf('FRAGMENT seed%d elapsed%.1f dv%.9f shots%d gated%d S2%d\n',cfg.seed,out.search_elapsed_s,r.total_dv_km_s,st.shots,st.gated,st.fragments);
 end
 if isfile(fullfile(folder,'STOP')), break; end
end
persist();
verifyStart=tic; out.independent=ctocscreen.propagateSchedule(s,p,true);
out.verification_elapsed_s=toc(verifyStart);
save(fullfile(folder,'result.mat'),'out','p','-v7.3');
if out.independent.passed
 elite=out; save(fullfile(folder,'elite.mat'),'elite','p','-v7.3');
end
fprintf('FRAGMENT_DONE seed%d search%.2f dv%.12f verified%d S2%d\n', ...
 cfg.seed,out.search_elapsed_s,r.total_dv_km_s,out.independent.passed,out.statistics.fragments);
 function persist()
  out.search_elapsed_s=elapsed+toc(start); out.cpu_seconds=cputime-cpuStart;
  out.random_state=stream.State;
  tmp=fullfile(folder,'checkpoint_next.mat'); save(tmp,'out','-v7.3'); movefile(tmp,checkpoint,'f');
 end
end
