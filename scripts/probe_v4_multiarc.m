function report=probe_v4_multiarc(label)
%PROBE_V4_MULTIARC Bounded historical component diagnostic, never a search.
% Reuse production events/STM/Newton equations; count gates before filtering.
if nargin<1, label='minimal_01'; end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
folder=fullfile(sim,'runs/v4/development',['multiarc_probe_' label]);
assert(~isfolder(folder),'Use a fresh label.'); mkdir(folder); clock=tic; budget=180;
c=ctocscreen.v4.defaults(); sig=ctocscreen.v4.signature();
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
files={'runs/v4/development/crossing_t2_w16_w16_defaults.mat', ...
 'runs/v4/development/incumbent_breadth_1800_seed888_01/report.mat', ...
 'runs/v4/development/incumbent_teammate_warm_1200_seed888_01/report.mat'};
hashes=cellfun(@(f)ctocscreen.v3FileHash(fullfile(sim,f)),files,'UniformOutput',false);
samples={}; old=load(fullfile(sim,files{1})); Q=old.report.q;
% Known positive control, selected before the probe from the archived table.
q=Q; q.tau=Q.tau(1:2); q.u=Q.u(1:2,:); q.T=Q.tau(3); q.witness(q.witness>q.T)=NaN;
p=Q; p.tau=Q.tau(1); p.u=Q.u(1,:); p.T=Q.tau(2); p.witness(p.witness>p.T)=NaN;
samples{1}=struct('source',files{1},'id',2,'kind','known_positive_control','q',q,'parent',p);
levels=[12 20 28 32];
for f=2:3
 z=load(fullfile(sim,files{f})); r=z.report; a=r.tree_audit;
 for level=levels
  candidates=find(a(:,2)>0&a(:,11)==0&a(:,4)==level);
  chosen=[];
  for candidateId=candidates(:).'
   par=a(candidateId,2);
   if a(candidateId,3)==a(par,3)+1&&a(candidateId,4)==a(par,4)+1, chosen=candidateId; break; end
  end
  assert(~isempty(chosen),'Missing predeclared sample stratum.');
  samples{end+1}=struct('source',files{f},'id',chosen,'kind','first_new_one_burn_one_visit', ...
   'q',r.node_controls{chosen},'parent',r.node_controls{a(chosen,2)}); %#ok<AGROW>
 end
end
manifest=struct('purpose','historical_component_diagnostic_not_cold_or_full_search','budget_s',budget, ...
 'signature',sig,'files',{files},'hashes',{hashes},'config',c,'levels',levels, ...
 'sampling','one known positive, then earliest eligible edge at each level in each run', ...
 'nonlinear','production Newton limits; one legacy and up to two missed candidates per node; boundary pair if crossed', ...
 'matlab',version);
save(fullfile(folder,'manifest.mat'),'manifest','samples'); records={};
for k=1:numel(samples)
 if toc(clock)>budget-20, break; end
 item=samples{k}; [a,tr]=ctocscreen.v4.replay(item.q,eph,c,false);
 assert(a.initial_passed&&a.height_passed&&~strcmp(a.status,'propagation_failure'));
 n=struct('q',tr.q,'actual',a,'trace',tr,'origin','diagnostic'); n.q.witness=a.witness_times_s;
 nr=inspect(n,item.parent); nr.sample=item; records{end+1}=nr; %#ok<AGROW>
 fprintf('SAMPLE %d id=%d visits=%d events=%d time_ok=%d dist_ok=%d cheap_missed=%d legacy=%s trials=%d\n', ...
  k,item.id,a.visit_count,nr.counts.events,nr.counts.time_ok,nr.counts.distance_ok, ...
  nr.counts.cheap_missed,nr.baseline.status,numel(nr.trials));
 report=struct('manifest',manifest,'records',{records},'elapsed_s',toc(clock));
 save(fullfile(folder,'report.mat'),'report');
end
report.samples_completed=numel(records); report.samples_planned=numel(samples);
report.source_unchanged=isequal(sig,ctocscreen.v4.signature());
report.inputs_unchanged=isequal(hashes,cellfun(@(f)ctocscreen.v3FileHash(fullfile(sim,f)),files,'UniformOutput',false));
assert(report.source_unchanged&&report.inputs_unchanged);
report.elapsed_s=toc(clock); save(fullfile(folder,'report.mat'),'report');
fprintf('DONE samples=%d/%d seconds=%.3f source_unchanged=%d\n',numel(records),numel(samples),toc(clock),report.source_unchanged);

 function out=inspect(node,parent)
  q=node.q; M=numel(q.tau); tau0=q.tau(end); u0=q.u(end,:).'; m=eph.model;
  [ev,ei]=ctocscreen.v4.events(node,eph,c); latest=m.horizon_s-max(0,sum(node.actual.distance_km>1)-2)*c.min_leg_s;
  tm=[ev.time]<=latest; dm=[ev.miss_km]<=c.absorb_miss_km;
  out=struct('counts',struct('events',numel(ev),'time_ok',sum(tm),'distance_ok',sum(dm), ...
   'both_ok',sum(tm&dm),'cheap_missed',0),'events_info',rmfield(ei,'sol'),'rows',[], ...
   'row_columns',{{'event','target','time','miss','time_ok','distance_ok','legacy_selected','extra','du_norm', ...
   'dtau','dtheta_j','dtheta_n','linear_residual','rcond_scaled','time_order_ok','parent_boundary_ok','scaled_step'}}, ...
   'baseline',[],'trials',{{}});
  [bc,br]=ctocscreen.v4.absorb(node,eph,c,min(5,max(.01,budget-20-toc(clock))));
  out.baseline=br;
  if ~isempty(bc), out.baseline.audit=audit(bc.q,node,parent); end
  if isempty(ev)||toc(clock)>budget-20, return; end
  valid=find(node.actual.distance_km<=1&isfinite(node.actual.witness_times_s));
  j=find(node.actual.distance_km<=1&abs(node.actual.witness_times_s-q.T)<1e-6,1);
  if isempty(j), [~,loc]=max(node.actual.witness_times_s(valid)); j=valid(loc); end
  tj=node.actual.witness_times_s(j); out.anchor_target=j;
  out.anchor_time=tj; out.parent_end=parent.T; out.original_departure=tau0;
  oldIds=find(tm&dm); [~,ord]=sort([ev(oldIds).miss_km]); oldIds=oldIds(ord);
  selected=[]; targets=[];
  for e=oldIds
   if ~ismember(ev(e).id,targets), selected(end+1)=e; targets(end+1)=ev(e).id; end %#ok<AGROW>
   if numel(selected)==c.absorb_candidates, break; end
  end
  eligible=find(tm); if isempty(eligible), return; end
  xp=ctocscreen.v4.stateAt(node.trace,tau0,'post');
  [~,~,sol]=ctocscreen.v3Arc(xp,tau0,max([ev(eligible).time]),m,c,true);
  rows=nan(numel(eligible),17); scale=[.1;.1;.1;3600;3600;3600];
  for kk=1:numel(eligible)
   e=eligible(kk); event=ev(e); [A,F]=equations(sol,j,tj,event.id,event.time,u0);
   dx=-((A.*scale.')\F).*scale;
   lo=0; if M>1, lo=q.tau(M-1); end
   t=[tau0+dx(4),tj+dx(5),event.time+dx(6)];
   order=t(1)>=lo&&t(2)>t(1)&&t(3)>=t(2)&&t(3)<=m.horizon_s;
   rows(kk,:)=[e,event.id,event.time,event.miss_km,true,dm(e),ismember(e,selected), ...
    norm(u0+dx(1:3))-norm(u0),norm(dx(1:3)),dx(4:6).',norm(A*dx+F), ...
    rcond(A.*scale.'),order,t(1)>=parent.T,max(abs(dx./scale))];
  end
  out.rows=rows;
  cheap=all(isfinite(rows(:,8:14)),2)&rows(:,8)<=c.absorb_max_extra_km_s&rows(:,13)<=c.absorb_residual_km;
  missed=cheap&rows(:,6)==0; out.counts.cheap_missed=sum(missed);
  % One lowest predicted extra cost and one smallest scaled correction.
  ids=find(missed); picks=[];
  if ~isempty(ids)
   [~,ix]=min(rows(ids,8)); picks=ids(ix);
   [~,ix]=min(rows(ids,17)); picks=unique([picks,ids(ix)],'stable');
  end
  for kk=picks
   if toc(clock)>budget-20, break; end
   row=rows(kk,:); e=ev(row(1));
   [qn,rr]=refine(e,false); trial=struct('row',row,'mode','legacy_time_bound','newton',rr,'audit',[]);
   if ~isempty(qn), trial.audit=audit(qn,node,parent); end
   out.trials{end+1}=trial;
   % A paired solve is needed only when the free solve actually crosses B.
   if rr.min_tau<parent.T-1e-6&&toc(clock)<budget-20
    [qn,rr]=refine(e,true); trial=struct('row',row,'mode','parent_time_bound','newton',rr,'audit',[]);
    if ~isempty(qn), trial.audit=audit(qn,node,parent); end
    out.trials{end+1}=trial;
   end
  end

  function [qn,info]=refine(e,safe)
   trialClock=tic; qn=[]; u=u0; tau=tau0; tJ=tj; tN=e.time;
   lo=0; if M>1, lo=q.tau(M-1); end
   if safe, lo=max(lo,parent.T); end
   xpre=ctocscreen.v4.stateAt(node.trace,tau0,'pre');
   info=struct('status','iteration_limit','iterations',0,'miss_km',[NaN NaN],'min_tau',tau0,'seconds',0);
   for it=1:c.absorb_newton_iterations
    info.iterations=it;
    if toc(trialClock)>5||toc(clock)>budget-15, info.status='budget'; break; end
    try
     if tau<=tau0, xm=ctocscreen.v4.stateAt(node.trace,tau,'pre'); else, xm=ctocscreen.v3Arc(xpre,tau0,tau,m,c); end
     [~,~,s2]=ctocscreen.v3Arc(xm+[0;0;0;u],tau,max(tJ,tN),m,c,true);
     [A,F]=equations(s2,j,tJ,e.id,tN,u);
    catch err
     info.status='propagation'; info.failure=err.identifier; break
    end
    info.miss_km=[norm(F(1:3)) norm(F(4:6))];
    if max(info.miss_km)<=c.absorb_newton_km
     qn=q; qn.u(M,:)=u.'; qn.tau(M)=tau; qn.T=tN; qn.witness(j)=tJ; qn.witness(e.id)=tN;
     info.status='converged'; break
    end
    sc=[.1;.1;.1;3600;3600;3600]; dx=-((A.*sc.')\F).*sc;
    dx=dx*min(1,min(.2/max(norm(dx(1:3)),eps),1800/max(abs(dx(4:6)))));
    u=u+dx(1:3); tau=min(max(tau+dx(4),lo),tJ-1); tJ=max(tau+1,tJ+dx(5)); tN=min(m.horizon_s,max(tJ,tN+dx(6)));
    info.min_tau=min(info.min_tau,tau);
   end
   info.seconds=toc(trialClock); info.final_tau=tau;
  end
 end
 function [A,F]=equations(sol,j,tj,n,tn,u)
  [rj,vj]=ctocscreen.v3QueryTargets(eph,j,tj); [rn,vn]=ctocscreen.v3QueryTargets(eph,n,tn);
  yj=deval(sol,tj); yn=deval(sol,tn); Pj=reshape(yj(7:end),6,6); Pn=reshape(yn(7:end),6,6);
  F=[yj(1:3)-rj.';yn(1:3)-rn.'];
  A=[Pj(1:3,4:6),-Pj(1:3,1:3)*u,yj(4:6)-vj.',zeros(3,1); ...
     Pn(1:3,4:6),-Pn(1:3,1:3)*u,zeros(3,1),yn(4:6)-vn.'];
 end
 function a=audit(q,node,parent)
  [v,~]=ctocscreen.v4.replay(q,eph,c,true); active=node.actual.distance_km<=1;
  added=find(v.distance_km<=1&~active);
  prefixBurns=numel(parent.tau); unchanged=numel(q.tau)>=prefixBurns&&isequal(q.x0,parent.x0) ...
   &&isequal(q.tau(1:prefixBurns),parent.tau)&&isequal(q.u(1:prefixBurns,:),parent.u);
  newTimes=q.tau(prefixBurns+1:end); prefixValid=unchanged&&all(newTimes>=parent.T);
  a=struct('q',q,'verification',v,'active_preserved',all(v.distance_km(active)<=1), ...
   'added',added,'extra',v.total_dv_km_s-node.actual.total_dv_km_s,'control_prefix_equal',unchanged, ...
   'parent_prefix_valid',prefixValid,'tau_shift',q.tau(end)-node.q.tau(end), ...
   'local_valid',v.independent&&v.initial_passed&&v.height_passed&&all(v.distance_km(active)<=1)&&~isempty(added));
 end
end
