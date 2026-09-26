function report=run_v4_bound4_pilot(label)
%RUN_V4_BOUND4_PILOT Cold, cost-bounded partial-coverage diagnostic (no rollout).
% Existing A-layer expansion; each physical node expanded once. Cached sibling
% queues and deferred paths support four-impulse backtracking without replaying
% previously tried edges. This is finite sampling, not complete beam-stack.
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
clock=tic; budget=180; reserve=18; depth=4; seed=888;
folder=fullfile(sim,'runs/v4/development',['bound4_' label]);
assert(~isfolder(folder),'Use a new output label.'); mkdir(folder);
c=ctocscreen.v4.defaults(struct('seed',seed,'budget_s',budget));
signature=ctocscreen.v4.signature();
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
stream=RandStream('mt19937ar','Seed',seed); memory=[];
manifest=struct('cold_start',true,'history_inputs',{{}},'initial_candidates',0,'resume_file','', ...
 'config',c,'signature',signature,'target_signature',eph.signature,'matlab',version, ...
 'total_budget_s',budget,'verification_reserve_s',reserve,'backtrack_burns',depth, ...
 'rollout',false,'joint_B',false,'absorb',false,'initial_root_count',4, ...
 'started_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss')));
save(fullfile(folder,'manifest.mat'),'manifest');
nodes={}; slots=repmat(struct('path',[],'deferred',{{}}),1,4);
seen=containers.Map('KeyType','char','ValueType','logical');
edgeTaken=false(1,0); events={}; roots={}; bestIds=[]; improvements=[];
stats=struct('expanded',0,'generated',0,'over_cost',0,'duplicate_children',0, ...
 'guidance_failures',0,'backtracks',0,'four_burn_backtracks',0,'short_backtracks',0, ...
 'deferred_resumes',0,'roots',0,'unique_edges',0,'expansion_seconds',0);
fid=fopen(fullfile(folder,'events.jsonl'),'w','n','UTF-8'); assert(fid>0);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
tick=0; deadline=budget-reserve;
while toc(clock)<deadline-0.1
 tick=tick+1; slot=1+mod(tick-1,4);
 if isempty(slots(slot).path)
  stats.roots=stats.roots+1;
  node=ctocscreen.v4.root(stats.roots,eph,c,stream);
  id=addNode(node,0); slots(slot).path=id;
  roots{end+1}=struct('id',node.root_id,'kind',node.root_kind,'x0',node.q.x0, ...
   'seed_target',node.seed_target,'created_s',toc(clock)); %#ok<AGROW>
  record('root',roots{end}); consider(id);
 end
 path=slots(slot).path; id=path(end); rec=nodes{id};
 if ~rec.expanded
  assert(rec.node.actual.total_dv_km_s<=c.search_max_dv_km_s,'Over-bound expansion.');
  available=deadline-toc(clock); if available<.1, break; end
  bt=tic;
  [kids,er,memory]=ctocscreen.v4.expand(rec.node,eph,c,stream,memory,min(c.action_seconds,available));
  stats.expansion_seconds=stats.expansion_seconds+toc(bt);
  stats.expanded=stats.expanded+1; stats.generated=stats.generated+numel(kids);
  stats.guidance_failures=stats.guidance_failures+er.failed_guidance;
  rec.expanded=true; scores=[]; childIds=[];
  for j=1:numel(kids)
   child=kids{j}; key=ctocscreen.v4.controlKey(child.q);
   if isKey(seen,key), stats.duplicate_children=stats.duplicate_children+1; continue; end
   seen(key)=true;
   assertPrefix(rec.node.q,child.q);
   if child.actual.total_dv_km_s>c.search_max_dv_km_s
    stats.over_cost=stats.over_cost+1; rec.over_cost=rec.over_cost+1;
    record('cost_rejected',struct('parent',id,'root',child.root_id,'visits',child.actual.visit_count, ...
     'dv',child.actual.total_dv_km_s,'day',child.q.T/86400));
    continue
   end
   if ~child.actual.initial_passed||~child.actual.height_passed ...
     ||strcmp(child.actual.status,'propagation_failure')||child.zero_gain>c.max_zero_gain
    continue
   end
   cid=addNode(child,id); consider(cid); childIds(end+1)=cid; %#ok<AGROW>
   scores(end+1,:)=[-child.actual.visit_count,child.actual.total_dv_km_s ...
    +child.actual.inclination_penalty+c.time_weight*35*max(0,child.q.T/864000-child.actual.visit_count/35)]; %#ok<AGROW>
  end
  if ~isempty(childIds), [~,ord]=sortrows(scores,[1 2]); rec.children=childIds(ord); end
  nodes{id}=rec;
  record('expand',struct('node',id,'root',rec.node.root_id,'visits',rec.node.actual.visit_count, ...
   'dv',rec.node.actual.total_dv_km_s,'admitted',numel(rec.children), ...
   'over_cost',rec.over_cost,'guidance_failures',er.failed_guidance,'seconds',er.seconds));
 end
 if rec.next<=numel(rec.children)
  cid=rec.children(rec.next); rec.next=rec.next+1; nodes{id}=rec;
  assert(~edgeTaken(cid),'Repeated branch traversal.'); edgeTaken(cid)=true;
  stats.unique_edges=stats.unique_edges+1;
  slots(slot).path=[path,cid];
 else
  rollback(slot);
 end
end
stats.search_end_s=toc(clock);
% Independent check of the best screened prefix, with fallback only if needed.
checked={}; winner=[];
for j=1:numel(bestIds)
 if toc(clock)>=budget-1, break; end
 candidate=nodes{bestIds(j)}.node; started=toc(clock);
 [v,~]=ctocscreen.v4.replay(candidate.q,eph,c,true);
 valid=v.independent&&v.initial_passed&&v.height_passed ...
  &&~strcmp(v.status,'propagation_failure')&&v.total_dv_km_s<=c.search_max_dv_km_s;
 entry=struct('node_id',bestIds(j),'q',candidate.q,'screened',candidate.actual, ...
  'verification',v,'partial_constraints_passed',valid,'started_s',started,'finished_s',toc(clock));
 checked{end+1}=entry; %#ok<AGROW>
 if valid&&(isempty(winner)||v.visit_count>winner.verification.visit_count ...
   ||(v.visit_count==winner.verification.visit_count&&v.total_dv_km_s<winner.verification.total_dv_km_s))
  winner=entry;
 end
 if valid&&v.visit_count==candidate.actual.visit_count, break; end
 % Preserve enough time for another check based on the observed check cost.
 if budget-toc(clock)<toc(clock)-started+1, break; end
end
stats.verification_seconds=toc(clock)-stats.search_end_s;
stats.source_unchanged=isequal(signature,ctocscreen.v4.signature());
assert(stats.source_unchanged);
stats.expansion_once_checked=stats.expanded==sum(cellfun(@(n)n.expanded,nodes));
stats.unique_edge_checked=stats.unique_edges==sum(edgeTaken);
assert(stats.expansion_once_checked&&stats.unique_edge_checked);
% Save compact tree audit; integration traces stay out of the archive.
audit=zeros(numel(nodes),9);
for j=1:numel(nodes)
 n=nodes{j}; audit(j,:)=[j,n.parent,n.node.root_id,numel(n.node.q.tau), ...
  n.node.actual.visit_count,n.node.actual.total_dv_km_s,n.node.q.T,n.expanded,n.next-1];
end
bestScreened=nodes{bestIds(1)}.node;
bestScreened=rmfield(bestScreened,'trace');
report=struct('manifest',manifest,'stats',stats,'roots',{roots},'tree_audit',audit, ...
 'tree_columns',{{'id','parent','root','burns','visits','dv','T','expanded','taken_children'}}, ...
 'events',{events},'improvements',improvements,'best_screened',bestScreened, ...
 'checked',{checked},'best',winner,'rng_state',stream.State,'finished',true);
report.stats.total_seconds=toc(clock); report.stats.within_budget=report.stats.total_seconds<=budget;
save(fullfile(folder,'report.mat'),'report');
report.stats.total_seconds=toc(clock); report.stats.within_budget=report.stats.total_seconds<=budget;
save(fullfile(folder,'report.mat'),'report');
if isempty(winner)
 fprintf('FINAL no independently checked admissible prefix total=%.3fs folder=%s\n',report.stats.total_seconds,folder);
else
 v=winner.verification;
 fprintf('FINAL independent=%d/35 J=%.12f partial_constraints=%d complete=%d T=%.6fd height=%.6f total=%.3fs within180=%d\n', ...
  v.visit_count,v.total_dv_km_s,winner.partial_constraints_passed,v.passed,winner.q.T/86400, ...
  v.min_altitude_lower_km,report.stats.total_seconds,report.stats.within_budget);
end
fprintf('STATS roots=%d expanded=%d overcost=%d backtracks=%d four=%d duplicate_children=%d edges=%d\n', ...
 stats.roots,stats.expanded,stats.over_cost,stats.backtracks,stats.four_burn_backtracks,stats.duplicate_children,stats.unique_edges);

 function id=addNode(node,parent)
  id=numel(nodes)+1;
  nodes{id}=struct('node',node,'parent',parent,'expanded',false,'children',[],'next',1,'over_cost',0);
  edgeTaken(id)=false; seen(ctocscreen.v4.controlKey(node.q))=true;
 end
 function consider(id)
  prior=[]; if ~isempty(bestIds), prior=bestIds(1); end
  pool=[bestIds,id]; rank=zeros(numel(pool),2);
  for k=1:numel(pool)
   a=nodes{pool(k)}.node.actual; rank(k,:)=[-a.visit_count,a.total_dv_km_s];
  end
  [~,ord]=sortrows(rank,[1 2]); bestIds=pool(ord(1:min(3,numel(ord))));
  if isempty(prior)||bestIds(1)~=prior
   a=nodes{bestIds(1)}.node.actual;
   improvements(end+1,:)=[toc(clock),a.visit_count,a.total_dv_km_s,nodes{bestIds(1)}.node.q.T,id];
   fprintf('BEST t=%.1f visits=%d J=%.9f root=%d day=%.5f\n',toc(clock),a.visit_count,a.total_dv_km_s, ...
    nodes{bestIds(1)}.node.root_id,nodes{bestIds(1)}.node.q.T/86400);
  end
 end
 function rollback(slot)
  path=slots(slot).path; current=nodes{path(end)};
  if numel(path)==1
   slots(slot).path=[];
   while ~isempty(slots(slot).deferred)
    saved=slots(slot).deferred{end}; slots(slot).deferred(end)=[];
    frame=nodes{saved(end)};
    if frame.next<=numel(frame.children)
     slots(slot).path=saved; stats.deferred_resumes=stats.deferred_resumes+1;
     record('resume_deferred',struct('node',saved(end),'root',frame.node.root_id)); return
    end
   end
   return
  end
  M=numel(current.node.q.tau); counts=cellfun(@(n)numel(n.node.q.tau),nodes(path));
  if M<depth, pos=1; else, pos=find(counts(1:end-1)<=M-depth,1,'last'); end
  assert(~isempty(pos));
  ancestor=nodes{path(pos)}; removed=M-numel(ancestor.node.q.tau);
  assert(removed==min(depth,M),'Rollback did not remove the requested burn count.');
  assertPrefix(ancestor.node.q,current.node.q);
  % Keep unfinished branches in the skipped frames for later service.
  for k=pos+1:numel(path)-1
   frame=nodes{path(k)};
   if frame.next<=numel(frame.children), slots(slot).deferred{end+1}=path(1:k); end
  end
  slots(slot).path=path(1:pos);
  stats.backtracks=stats.backtracks+1;
  if removed==depth, stats.four_burn_backtracks=stats.four_burn_backtracks+1;
  else, stats.short_backtracks=stats.short_backtracks+1; end
  record('backtrack',struct('from',path(end),'to',path(pos),'root',current.node.root_id, ...
   'removed_burns',removed,'from_day',current.node.q.T/86400,'to_day',ancestor.node.q.T/86400, ...
   'before_visits',current.node.actual.visit_count,'after_visits',ancestor.node.actual.visit_count, ...
   'cost_rejections_here',current.over_cost));
 end
 function record(kind,payload)
  row=struct('kind',kind,'elapsed_s',toc(clock),'payload',payload);
  events{end+1}=row; fprintf(fid,'%s\n',jsonencode(row));
 end
end
function assertPrefix(p,q)
M=numel(p.tau);
assert(numel(q.tau)>=M&&isequal(p.x0,q.x0)&&isequal(p.tau,q.tau(1:M)) ...
 &&isequal(p.u,q.u(1:M,:)),'Physical prefix changed.');
end
