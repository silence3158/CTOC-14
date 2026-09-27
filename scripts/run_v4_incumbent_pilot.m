function report=run_v4_incumbent_pilot(label,budget,seed,options)
%RUN_V4_INCUMBENT_PILOT Cold two-phase test of incumbent-bounded backtracking.
% Phase 1: the formal V4 search (search.m) until its first independently
% verified complete solution (stop_on_complete, capped). Phase 2: depth-first
% branch-and-bound that starts from that incumbent's own control path; the
% bound is the best independently verified complete cost and tightens on
% every improvement (anytime DFBnB; Hansen & Zhou 2007, sec. 4.3). Rollback
% removes 4 actual impulses as in run_v4_bound4_pilot. The broaden_at-th actual
% landing on a node expands it again with its tried encounters excluded; the
% return_limit-th landing stops it and cascades (breadth cutoff, after
% Ginsberg & Harvey 1990). Proposals whose predicted two-body leg cost exceeds
% the remaining bound plus a margin are dropped before J2 guidance.
% Empty historical input: phase 2 only uses this run's phase-1 result.
if nargin<1||isempty(label), label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
if nargin<2||isempty(budget), budget=300; end
if nargin<3||isempty(seed), seed=888; end
if nargin<4||isempty(options), options=struct(); end
o=struct('phase1_cap_s',200,'broaden_at',2,'return_limit',4,'depth',4,'margin_km_s',0.1, ...
 'reserve_s',18,'epsilon_km_s',1e-6,'min_phase2_s',30,'exclude_window_s',600);
names=fieldnames(options); assert(isempty(names)||all(isfield(o,names)),'ctocscreen:v4:option','Unknown option.');
for kn=1:numel(names), o.(names{kn})=options.(names{kn}); end
validateattributes(budget,{'double'},{'scalar','finite','positive'});
validateattributes(seed,{'double'},{'scalar','finite','integer','>=',0,'<=',2^32-1});
assert(o.broaden_at>=1&&o.broaden_at<o.return_limit&&o.depth>=1&&o.margin_km_s>=0);
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
clock=tic;
folder=fullfile(sim,'runs/v4/development',['incumbent_' label]);
assert(~isfolder(folder),'Use a new output label.'); mkdir(folder);
c=ctocscreen.v4.defaults(struct('seed',seed,'budget_s',budget));
signature=ctocscreen.v4.signature();
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
horizon=eph.model.horizon_s; deadline=budget-o.reserve_s; eps0=o.epsilon_km_s;
experimentFiles={[mfilename('fullpath') '.m'],fullfile(sim,'scripts/v4_bound4_return_landings.m')};
experimentText=cellfun(@fileread,experimentFiles,'UniformOutput',false);
save(fullfile(folder,'experiment_source.mat'),'experimentFiles','experimentText');
manifest=struct('cold_start',true,'history_inputs',{{}},'initial_candidates',0,'resume_file','', ...
 'config',c,'options',o,'signature',signature,'target_signature',eph.signature,'matlab',version, ...
 'total_budget_s',budget,'seed',seed,'phase1','formal V4 search.m with stop_on_complete', ...
 'phase2','incumbent-bounded DFBnB from the incumbent path; broaden/limit landings; budget filter', ...
 'started_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss')));
save(fullfile(folder,'manifest.mat'),'manifest');
fid=fopen(fullfile(folder,'events.jsonl'),'w','n','UTF-8'); assert(fid>0);
cleanup=onCleanup(@()fclose(fid)); events={};
%% Phase 1: formal search until its first independently verified complete solution.
cap=min(o.phase1_cap_s,deadline-o.min_phase2_s-toc(clock)); assert(cap>=60,'Budget too small for phase 1.');
c1=ctocscreen.v4.defaults(struct('seed',seed,'budget_s',cap,'stop_on_complete',1));
t1=toc(clock); result1=ctocscreen.v4.search(fullfile(folder,'phase1_beam'),eph,c1);
phase1=struct('cap_s',cap,'started_s',t1,'finished_s',toc(clock), ...
 'found_complete',~isempty(result1.best_verified),'first_complete_search_s',result1.stats.first_complete_s, ...
 'J',NaN,'visits',result1.verification.visit_count,'stats',result1.stats);
if phase1.found_complete, phase1.J=result1.best_verified.actual.total_dv_km_s; end
record('phase1',rmfield(phase1,'stats'));
fprintf('PHASE1 complete=%d J=%.9f finished=%.1fs\n',phase1.found_complete,phase1.J,phase1.finished_s);
%% Phase 2 state.
stream=RandStream('mt19937ar','Seed',seed); stream.State=result1.rng_state; memory=[];
nodes={}; seen=containers.Map('KeyType','char','ValueType','logical');
edgeTaken=false(1,0); returnCounts=zeros(1,0); stopped=false(1,0);
path=[]; deferred={}; suppressedPaths={}; rootIndex=result1.stats.root_count;
rootRecords={}; snapshots={}; snapshotLimits=[300 600 1200]; snapshotNext=1;
stats=struct('expansions',0,'broadenings',0,'generated',0,'admitted',0,'broaden_admitted',0, ...
 'over_bound',0,'duplicates',0,'budget_filtered',0,'excluded',0,'guided',0,'guidance_failures',0, ...
 'bound_skipped',0,'backtracks',0,'four_burn_backtracks',0,'short_backtracks',0, ...
 'threshold_hits',0,'threshold_root_retirements',0,'exhausted_root_retirements',0, ...
 'deferred_resumes',0,'suppressed_deferred_paths',0,'verifications',0,'verify_failed',0, ...
 'improvements',0,'fresh_roots',0,'unique_edges',0,'virtual_edges',0,'expansion_seconds',0, ...
 'verification_seconds',0,'dead_rollbacks',0,'complete_rollbacks',0,'path_build_seconds',0);
incumbent=[]; bound=Inf; improvements=zeros(0,7); boundHistory=zeros(0,2); M1=NaN;
if phase1.found_complete
 source=result1.best_verified; bound=source.actual.total_dv_km_s; M1=numel(source.q.tau);
 incumbent=struct('q',source.q,'verification',source.actual,'node',NaN,'source','phase1_beam', ...
  'found_s',phase1.finished_s,'from_broaden',false,'revised_from_end',NaN);
else
 source=result1.best; % fallback: deepest partial, unbounded until a complete appears
end
boundHistory(end+1,:)=[toc(clock),bound];
pb=tic; path=buildPath(source); stats.path_build_seconds=toc(pb);
rememberRoot(nodes{path(1)}.node,path(1),'phase1_lineage');
record('phase2_start',struct('bound',bound,'path_nodes',numel(path),'burns',numel(source.q.tau)));
%% Phase 2: incumbent-bounded backtracking (single lineage; fresh roots only after it retires).
while toc(clock)<deadline-0.1
 if budget>300&&snapshotNext<=numel(snapshotLimits)&&toc(clock)>=snapshotLimits(snapshotNext)
  saveSnapshot(snapshotLimits(snapshotNext),false); snapshotNext=snapshotNext+1;
 end
 if isempty(path)
  rootIndex=rootIndex+1; stats.fresh_roots=stats.fresh_roots+1;
  node=ctocscreen.v4.root(rootIndex,eph,c,stream); id=addNode(node,0,false); path=id;
  rememberRoot(node,id,'fresh_phase2');
  record('root',struct('node',id,'root',rootIndex,'kind',node.root_kind,'bound',bound));
 end
 id=path(end); rec=nodes{id}; assert(~any(stopped(path)),'Stopped subtree reentered.');
 a=rec.node.actual;
 if a.passed||rec.verify_failed, stats.complete_rollbacks=stats.complete_rollbacks+1; rollback(); continue; end
 if a.total_dv_km_s>=bound-eps0, stats.dead_rollbacks=stats.dead_rollbacks+1; rollback(); continue; end
 if rec.node.q.T>=horizon-1, rollback(); continue; end
 if ~rec.expanded
  if ~expandNode(id,false), break; end
 elseif rec.broaden_pending
  if ~expandNode(id,true), break; end
 end
 rec=nodes{id}; taken=false;
 while rec.next<=numel(rec.children)
  cid=rec.children(rec.next); rec.next=rec.next+1;
  if nodes{cid}.node.actual.total_dv_km_s>=bound-eps0, stats.bound_skipped=stats.bound_skipped+1; continue; end
  assert(~edgeTaken(cid),'Repeated branch traversal.'); edgeTaken(cid)=true;
  stats.unique_edges=stats.unique_edges+1; path=[path,cid]; taken=true; break %#ok<AGROW>
 end
 nodes{id}=rec;
 if ~taken, rollback(); end
end
stats.phase2_end_s=toc(clock);
if budget>300, saveSnapshot(budget,true); end
%% Final result: the incumbent carries its own independent verification.
if isempty(incumbent)
 [~,jb]=max(cellfun(@(n)n.node.actual.visit_count-1e-3*n.node.actual.total_dv_km_s,nodes));
 [v,~]=ctocscreen.v4.replay(nodes{jb}.node.q,eph,c,true);
 final=struct('q',nodes{jb}.node.q,'verification',v,'node',jb,'source','partial_no_complete');
else
 final=incumbent;
end
stats.source_unchanged=isequal(signature,ctocscreen.v4.signature());
stats.experiment_source_unchanged=isequal(experimentText,cellfun(@fileread,experimentFiles,'UniformOutput',false));
assert(stats.source_unchanged&&stats.experiment_source_unchanged,'Sources changed during this experiment.');
% Audits recomputed from the event log.
landingAudit=zeros(size(returnCounts)); broadenLandings=[];
for ie=1:numel(events)
 if strcmp(events{ie}.kind,'backtrack'), lt=events{ie}.payload.to; landingAudit(lt)=landingAudit(lt)+1; end
 if strcmp(events{ie}.kind,'expand')&&events{ie}.payload.broaden, broadenLandings(end+1)=events{ie}.payload.landings; end %#ok<AGROW>
end
stats.return_count_checked=isequal(landingAudit,returnCounts)&&all(returnCounts<=o.return_limit);
stats.broaden_landing_checked=all(broadenLandings==o.broaden_at);
expCount=sum(cellfun(@(n)double(n.expanded)+n.broadened,nodes));
stats.expansion_count_checked=expCount==stats.expansions;
stats.unique_edge_checked=sum(edgeTaken)==stats.unique_edges+stats.virtual_edges;
assert(stats.return_count_checked&&stats.broaden_landing_checked&&stats.expansion_count_checked&&stats.unique_edge_checked);
audit=zeros(numel(nodes),15); controls=cell(1,numel(nodes)); structures=cell(1,numel(nodes));
for ia=1:numel(nodes)
 na=nodes{ia}; audit(ia,:)=[ia,na.parent,numel(na.node.q.tau),na.node.actual.visit_count,na.node.actual.total_dv_km_s, ...
  na.node.q.T,na.expanded,na.broadened,numel(na.children),na.next-1,na.virtual,returnCounts(ia), ...
  na.node.root_id,na.common_prefix,na.created_s];
 controls{ia}=na.node.q; structures{ia}=visitStructure(na.node);
end
report=struct('manifest',manifest,'options',o,'phase1',phase1,'stats',stats,'final',final, ...
 'phase1_incumbent_J',phase1.J,'improvements',improvements, ...
 'improvement_columns',{{'t_s','J','node','from_broaden','root','revised_from_end','burns'}}, ...
 'bound_history',boundHistory,'tree_audit',audit, ...
 'tree_columns',{{'id','parent','burns','visits','dv','T','expanded','broadened','children','taken','virtual','landings','root','common_prefix','created_s'}}, ...
 'root_records',{rootRecords},'snapshots',{snapshots},'node_controls',{controls},'visit_structures',{structures}, ...
 'events',{events},'stopped_nodes',find(stopped),'suppressed_paths',{suppressedPaths},'rng_state',stream.State);
report.stats.total_seconds=toc(clock); report.stats.within_budget=report.stats.total_seconds<=budget;
save(fullfile(folder,'report.mat'),'report');
v=final.verification;
if v.passed
 verified=final; save(fullfile(folder,'verified_complete.mat'),'verified');
 if v.total_dv_km_s<=c.search_max_dv_km_s, save(fullfile(folder,'elite.mat'),'verified');
 else, save(fullfile(folder,'rejected_complete.mat'),'verified'); end
end
report.stats.total_seconds=toc(clock); report.stats.within_budget=report.stats.total_seconds<=budget;
save(fullfile(folder,'report.mat'),'report');
fprintf('FINAL independent=%d/35 J=%.12f passed=%d phase1_J=%.9f improvements=%d total=%.3fs withinBudget=%d\n', ...
 v.visit_count,v.total_dv_km_s,v.passed,phase1.J,stats.improvements,report.stats.total_seconds,report.stats.within_budget);
fprintf('STATS expansions=%d broadenings=%d admitted=%d (%.2f/exp) over_bound=%d filtered=%d excluded=%d backtracks=%d limit_hits=%d\n', ...
 stats.expansions,stats.broadenings,stats.admitted,stats.admitted/max(1,stats.expansions),stats.over_bound, ...
 stats.budget_filtered,stats.excluded,stats.backtracks,stats.threshold_hits);

 function ids=buildPath(src)
  % Path nodes P_k: initial state and the first k impulses, ending at the last
  % visit of the arc after impulse k (else just before impulse k+1; the root
  % ends at epoch). Exact prefix reuse keeps the replay identical.
  q=src.q; w=src.actual.witness_times_s; vis=src.actual.distance_km<=1; M=numel(q.tau);
  edges=[0;q.tau(:);q.T]; ids=zeros(1,M+1); prev=[]; prevId=0;
  for kp=0:M
   qk=q; qk.tau=q.tau(1:kp); qk.u=q.u(1:kp,:);
   if kp<M
    wk=w(vis&w>edges(kp+1)&w<=edges(kp+2));
    if ~isempty(wk), Tk=max(wk); elseif kp==0, Tk=0; else, Tk=edges(kp+2); end
   else
    Tk=q.T;
   end
   qk.T=Tk; qk.witness=nan(35,1); hint=vis&w<=Tk; qk.witness(hint)=w(hint);
   [ak,tk]=ctocscreen.v4.replay(qk,eph,c,false,prev);
   assert(~strcmp(ak.status,'propagation_failure')&&ak.initial_passed&&ak.height_passed,'Incumbent prefix replay failed.');
   node=src; node.q=tk.q; node.q.witness=ak.witness_times_s; node.actual=ak; node.trace=tk;
   node.attempts=0; node.zero_gain=0; node.generation=kp; node.origin='incumbent_prefix'; node.heuristic_H=NaN;
   pid=addNode(node,prevId,false); nodes{pid}.virtual=true;
   if prevId>0
    nodes{prevId}.children=pid; nodes{prevId}.next=2; edgeTaken(pid)=true; stats.virtual_edges=stats.virtual_edges+1;
    newly=find(ak.distance_km<=1&prev.actual.distance_km>1);
    nodes{prevId}.tried=[newly(:),ak.witness_times_s(newly)];
   end
   ids(kp+1)=pid; prev=node; prevId=pid;
  end
  if isfinite(bound), assert(nodes{ids(end)}.node.actual.passed,'Incumbent replay is not complete.'); end
 end
 function ok=expandNode(id,broaden)
  ok=false; rec=nodes{id}; available=deadline-toc(clock);
  if available<0.5, return; end
  ok=true; eo=struct('max_leg_dv',Inf,'exclude',rec.tried,'exclude_window_s',o.exclude_window_s);
  if isfinite(bound), eo.max_leg_dv=bound-rec.node.actual.total_dv_km_s+o.margin_km_s; end
  bt=tic;
  [kids,er,memory]=ctocscreen.v4.expand(rec.node,eph,c,stream,memory,min(c.action_seconds,available),eo);
  stats.expansion_seconds=stats.expansion_seconds+toc(bt); stats.expansions=stats.expansions+1;
  stats.generated=stats.generated+numel(kids); stats.budget_filtered=stats.budget_filtered+er.budget_filtered;
  stats.excluded=stats.excluded+er.excluded; stats.guided=stats.guided+er.guided;
  stats.guidance_failures=stats.guidance_failures+er.failed_guidance;
  rec.tried=[rec.tried;er.tried]; lineage=broaden||rec.from_broaden;
  childIds=zeros(1,0); scores=zeros(0,2); overBound=0; completes=0;
  for jk=1:numel(kids)
   kid=kids{jk}; key=ctocscreen.v4.controlKey(kid.q);
   if isKey(seen,key), stats.duplicates=stats.duplicates+1; continue; end
   seen(key)=true; assertPrefix(rec.node.q,kid.q); J=kid.actual.total_dv_km_s;
   if J>=bound-eps0, overBound=overBound+1; continue; end
   if ~kid.actual.initial_passed||~kid.actual.height_passed||strcmp(kid.actual.status,'propagation_failure') ...
     ||kid.zero_gain>c.max_zero_gain
    continue
   end
   cid=addNode(kid,id,lineage); childIds(end+1)=cid; %#ok<AGROW>
   % Scalar first: inside [ ], a continued line starting with " +x" would
   % become a separate element (the ordering bug found in run_v4_bound4_pilot).
   f=J+kid.actual.inclination_penalty+c.time_weight*35*max(0,kid.q.T/864000-kid.actual.visit_count/35);
   scores(end+1,:)=[-kid.actual.visit_count,f]; %#ok<AGROW>
   if kid.actual.passed, completes=completes+1; checkComplete(cid); end
  end
  if ~isempty(childIds), [~,ord]=sortrows(scores,[1 2]); rec.children=[rec.children,childIds(ord)]; end
  stats.over_bound=stats.over_bound+overBound; stats.admitted=stats.admitted+numel(childIds);
  if broaden
   rec.broaden_pending=false; rec.broadened=rec.broadened+1;
   stats.broadenings=stats.broadenings+1; stats.broaden_admitted=stats.broaden_admitted+numel(childIds);
  end
  landings=returnCounts(id); rec.expanded=true; nodes{id}=rec;
  record('expand',struct('node',id,'root',rec.node.root_id,'common_prefix',rec.common_prefix, ...
   'broaden',broaden,'landings',landings,'virtual',rec.virtual, ...
   'burns',numel(rec.node.q.tau),'visits',rec.node.actual.visit_count,'dv',rec.node.actual.total_dv_km_s, ...
   'bound',bound,'max_leg_dv',eo.max_leg_dv,'generated',numel(kids),'admitted',numel(childIds), ...
   'over_bound',overBound,'budget_filtered',er.budget_filtered,'excluded',er.excluded,'guided',er.guided, ...
   'guidance_failures',er.failed_guidance,'completes',completes,'seconds',er.seconds));
 end
 function checkComplete(cid)
  n=nodes{cid}.node; J=n.actual.total_dv_km_s;
  if J>=bound-eps0, return; end
  vt=tic; [v,~]=ctocscreen.v4.replay(n.q,eph,c,true);
  stats.verification_seconds=stats.verification_seconds+toc(vt); stats.verifications=stats.verifications+1;
  valid=v.passed&&v.independent&&v.initial_passed&&v.height_passed&&v.total_dv_km_s<bound-eps0;
  if ~valid
   nodes{cid}.verify_failed=true; stats.verify_failed=stats.verify_failed+1;
   record('verify_failed',struct('node',cid,'visits',v.visit_count,'dv',v.total_dv_km_s,'status',v.status)); return
  end
  revised=revisedFromEnd(cid);
  incumbent=struct('q',n.q,'verification',v,'node',cid,'source','phase2_backtrack','found_s',toc(clock), ...
   'from_broaden',nodes{cid}.from_broaden,'revised_from_end',revised);
  bound=v.total_dv_km_s; boundHistory(end+1,:)=[toc(clock),bound]; stats.improvements=stats.improvements+1;
  improvements(end+1,:)=[toc(clock),bound,cid,nodes{cid}.from_broaden,n.root_id,revised,numel(n.q.tau)];
  improved=incumbent; save(fullfile(folder,sprintf('improvement_%03d.mat',stats.improvements)),'improved');
  record('improvement',struct('node',cid,'dv',bound,'from_broaden',nodes{cid}.from_broaden, ...
   'revised_from_end',revised,'burns',numel(n.q.tau)));
  fprintf('IMPROVED t=%.1f J=%.9f node=%d revised_from_end=%g broaden_lineage=%d\n', ...
   toc(clock),bound,cid,revised,nodes{cid}.from_broaden);
  if bound<c.notify_dv_km_s
   fprintf('V4 USER_ACCEPTANCE_READY independent=35/35 raw_dv=%.12f time=%.3f\n',bound,toc(clock));
   witness=incumbent; save(fullfile(folder,'acceptance_ready.mat'),'witness','manifest');
  end
 end
 function r=revisedFromEnd(cid)
  % Impulses of the phase-1 incumbent replaced by this lineage (NaN for fresh roots).
  r=NaN; anc=cid;
  while anc>0&&~nodes{anc}.virtual, anc=nodes{anc}.parent; end
  if anc>0&&isfinite(M1), r=M1-numel(nodes{anc}.node.q.tau); end
 end
 function rollback()
  if isempty(path), return; end
  if isscalar(path)
   root=path; path=[];
   while ~isempty(deferred)
    saved=deferred{end}; deferred(end)=[];
    if any(stopped(saved)), suppress(saved); continue; end
    frame=nodes{saved(end)};
    if ~frame.expanded||frame.next<=numel(frame.children)
     path=saved; stats.deferred_resumes=stats.deferred_resumes+1;
     record('resume_deferred',struct('node',saved(end))); return
    end
   end
   stats.exhausted_root_retirements=stats.exhausted_root_retirements+1;
   record('root_retired',struct('node',root,'reason','finite_candidates_exhausted')); return
  end
  counts=cellfun(@(n)numel(n.node.q.tau),nodes(path)); oldPath=path;
  [path,returnCounts,stopped,legs,hits]=v4_bound4_return_landings(oldPath,counts,returnCounts,stopped,o.depth,o.return_limit);
  for leg=1:size(legs,1)
   from=legs(leg,1); to=legs(leg,2); removed=legs(leg,3);
   assertPrefix(nodes{to}.node.q,nodes{from}.node.q);
   pos=find(oldPath==to,1); fromPos=find(oldPath==from,1);
   for kd=pos+1:fromPos-1
    frame=nodes{oldPath(kd)};
    if ~frame.expanded||frame.next<=numel(frame.children), deferred{end+1}=oldPath(1:kd); end %#ok<AGROW>
   end
   stats.backtracks=stats.backtracks+1;
   if removed==o.depth, stats.four_burn_backtracks=stats.four_burn_backtracks+1;
   else, stats.short_backtracks=stats.short_backtracks+1; end
   record('backtrack',struct('from',from,'to',to,'removed_burns',removed,'return_count',legs(leg,4), ...
    'to_burns',numel(nodes{to}.node.q.tau),'to_dv',nodes{to}.node.actual.total_dv_km_s,'bound',bound,'cascade',leg>1));
   if ismember(to,hits)
    stats.threshold_hits=stats.threshold_hits+1;
    record('return_limit',struct('node',to,'returns',legs(leg,4),'is_root',nodes{to}.parent==0));
   end
  end
  if ~isempty(hits)
   keep=true(size(deferred));
   for ks=1:numel(keep), if any(stopped(deferred{ks})), keep(ks)=false; suppress(deferred{ks}); end, end
   deferred=deferred(keep);
  end
  if isempty(path)
   stats.threshold_root_retirements=stats.threshold_root_retirements+1;
   record('root_retired',struct('node',legs(end,2),'reason','return_limit')); return
  end
  dest=path(end);
  if returnCounts(dest)==o.broaden_at&&nodes{dest}.expanded&&nodes{dest}.node.actual.total_dv_km_s<bound-eps0
   nodes{dest}.broaden_pending=true; record('broaden_scheduled',struct('node',dest,'returns',returnCounts(dest)));
  end
 end
 function id=addNode(node,parent,fromBroaden)
  id=numel(nodes)+1;
  nodes{id}=struct('node',node,'parent',parent,'expanded',false,'children',zeros(1,0),'next',1, ...
   'broaden_pending',false,'broadened',0,'from_broaden',fromBroaden,'tried',zeros(0,2), ...
   'verify_failed',false,'virtual',false,'common_prefix',commonPrefix(node.q),'created_s',toc(clock));
  edgeTaken(id)=false; seen(ctocscreen.v4.controlKey(node.q))=true; returnCounts(id)=0; stopped(id)=false;
 end
 function n=commonPrefix(q)
  % Exact shared controls with this run's phase-1 source; -1 means new x0.
  n=-1; if ~isequal(q.x0,source.q.x0), return; end
  n=0;
  for j=1:min(numel(q.tau),numel(source.q.tau))
   if q.tau(j)~=source.q.tau(j)||~isequal(q.u(j,:),source.q.u(j,:)), break; end
   n=j;
  end
 end
 function s=visitStructure(node)
  ids=find(node.actual.distance_km<=1); times=node.actual.witness_times_s(ids);
  [times,order]=sort(times); ids=ids(order); arcs=zeros(size(ids));
  for j=1:numel(ids), arcs(j)=sum(node.q.tau<times(j)); end
  s=struct('ids',ids,'times',times,'arcs',arcs);
 end
 function rememberRoot(node,id,origin)
  entry=struct('node',id,'root',node.root_id,'kind',node.root_kind,'origin',origin, ...
   'x0',node.q.x0,'initial_orbit',node.actual.initial_orbit,'created_s',nodes{id}.created_s);
  rootRecords{end+1}=entry; record('root_detail',entry);
 end
 function saveSnapshot(limit,isFinal)
  counts=zeros(numel(rootRecords),7);
  for j=1:numel(rootRecords)
   ids=find(cellfun(@(n)n.node.root_id==rootRecords{j}.root,nodes));
   active=ids(cellfun(@(n)n.expanded,nodes(ids)));
   minPrefix=NaN; if ~isempty(active), minPrefix=min(cellfun(@(n)n.common_prefix,nodes(active))); end
   counts(j,:)=[rootRecords{j}.root,numel(ids),numel(active), ...
    sum(cellfun(@(n)n.broadened,nodes(ids))), ...
    max(cellfun(@(n)n.node.actual.visit_count,nodes(ids))),minPrefix,returnCounts(rootRecords{j}.node)];
  end
  snapshot=struct('limit_s',limit,'actual_s',toc(clock),'final',isFinal, ...
   'incumbent',incumbent,'stats',stats,'roots',{rootRecords},'root_summary',counts, ...
   'root_columns',{{'root','nodes','expanded_nodes','broadenings','max_screened_visits','min_shared_prefix','root_returns'}});
  snapshots{end+1}=snapshot;
  save(fullfile(folder,sprintf('snapshot_%04d.mat',limit)),'snapshot');
  record('snapshot',struct('limit_s',limit,'actual_s',snapshot.actual_s,'bound',bound, ...
   'fresh_roots',stats.fresh_roots,'expansions',stats.expansions,'final',isFinal));
  fprintf('SNAPSHOT limit=%d actual=%.2f J=%.9f fresh_roots=%d expansions=%d\n', ...
   limit,snapshot.actual_s,bound,stats.fresh_roots,stats.expansions);
 end
 function suppress(saved)
  suppressedPaths{end+1}=saved; stats.suppressed_deferred_paths=stats.suppressed_deferred_paths+1;
  record('suppress_deferred',struct('node',saved(end),'blocked_by',saved(stopped(saved))));
 end
 function record(kind,payload)
  row=struct('kind',kind,'elapsed_s',toc(clock),'payload',payload);
  events{end+1}=row; fprintf(fid,'%s\n',jsonencode(row));
 end
end
function assertPrefix(p,q)
M=numel(p.tau);
assert(numel(q.tau)>=M,'Physical prefix changed.');
assert(isequal(p.x0,q.x0)&&isequal(p.tau(:),reshape(q.tau(1:M),[],1))&&isequal(p.u,q.u(1:M,:)), ...
 'Physical prefix changed.');
end
