function report=eval_v4_backtracking(label,secondsPerDepth)
%EVAL_V4_BACKTRACKING Authorized historical diagnostic, NOT a cold search.
% Remove 1/3/6 final impulses from the SAME source trajectory; reconstruct
% every affected suffix from its real prefix state, using bounded DFS over
% existing V4 expansion proposals. No original suffix states are restored.
% This is a sampled diagnostic, not the complete Beam-Stack algorithm.
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
if nargin<2, secondsPerDepth=60; end
validateattributes(secondsPerDepth,{'double'},{'scalar','positive','finite'});
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
outdir=fullfile(sim,'runs/v4/development',['backtrack_' label]);
assert(~isfolder(outdir),'Output exists; use a new label.'); mkdir(outdir);
source=fullfile(sim,'runs/v4/search/v4_three_steps_hev_900_seed888/verified_complete.mat');
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(struct('branch_count',4));
s=load(source,'verified'); original=s.verified; q0=original.q;
signature=ctocscreen.v4.signature(); totalClock=tic;
[baseline,~]=ctocscreen.v4.replay(q0,eph,c,true);
assert(baseline.passed&&baseline.independent,'Source failed independent replay.');
report=struct('diagnostic_only',true,'historical_input_authorized',true,'cold_start',false, ...
 'source',source,'source_sha256',ctocscreen.v3FileHash(source),'config',c, ...
 'signature',signature,'matlab',version,'baseline_q',q0,'baseline',baseline, ...
 'depths',[1 3 6],'seconds_per_depth',secondsPerDepth,'groups',{{}},'finished',false);
fprintf('BASE independent=%d visits=%d J=%.12f M=%d T=%.6f days\n', ...
 baseline.passed,baseline.visit_count,baseline.total_dv_km_s,numel(q0.tau),q0.T/86400);
save(fullfile(outdir,'report.mat'),'report');
for depth=report.depths
 groupClock=tic; first=numel(q0.tau)-depth+1; p=q0; p.T=q0.tau(first);
 p.tau=q0.tau(1:first-1); p.u=q0.u(1:first-1,:); p.witness(:)=NaN;
 % Retain only independently rechecked witnesses at/before the cut.
 keep=baseline.witness_times_s<=p.T; p.witness(keep)=baseline.witness_times_s(keep);
 [pa,pt]=ctocscreen.v4.replay(p,eph,c);
 assert(pa.initial_passed&&pa.height_passed&&~strcmp(pa.status,'propagation_failure'));
 root=original; root.q=pt.q; root.q.witness=pa.witness_times_s;
 root.actual=pa; root.trace=pt; root.attempts=0; root.zero_gain=0;
 root.generation=0; root.heuristic_H=NaN; root.origin='authorized_historical_prefix';
 stream=RandStream('mt19937ar','Seed',927001); memory=[]; stack={root};
 seen=containers.Map('KeyType','char','ValueType','logical'); candidates={};
 trials=struct('elapsed_s',{},'visits',{},'dv',{},'mission_day',{},'new_burns',{},'complete',{});
 expansions=0; childrenCount=0; duplicates=0; boundPruned=0; restarts=0;
 guidanceFailures=0; completeCount=0; bestPartial=pa.visit_count; searchClock=tic;
 fprintf('CUT depth=%d day=%.6f retained=%d prefixJ=%.9f revoked=%s\n', ...
 depth,p.T/86400,pa.visit_count,pa.total_dv_km_s,mat2str(find(pa.distance_km>1).'));
 while toc(searchClock)<secondsPerDepth
  if isempty(stack), stack={root}; restarts=restarts+1; end
  parent=stack{end}; stack(end)=[];
  available=secondsPerDepth-toc(searchClock); if available<.1, break; end
  [kids,er,memory]=ctocscreen.v4.expand(parent,eph,c,stream,memory,min(c.action_seconds,available));
  expansions=expansions+1; childrenCount=childrenCount+numel(kids);
  guidanceFailures=guidanceFailures+er.failed_guidance; open={}; scores=[];
  for j=1:numel(kids)
   node=kids{j}; key=ctocscreen.v4.controlKey(node.q);
   if isKey(seen,key), duplicates=duplicates+1; continue; end
   seen(key)=true;
   assert(isequal(node.q.x0,p.x0)&&isequal(node.q.tau(1:first-1),p.tau) ...
    &&isequal(node.q.u(1:first-1,:),p.u),'Fixed prefix changed.');
   bestPartial=max(bestPartial,node.actual.visit_count);
   trials(end+1)=struct('elapsed_s',toc(searchClock),'visits',node.actual.visit_count, ...
    'dv',node.actual.total_dv_km_s,'mission_day',node.q.T/86400, ...
    'new_burns',numel(node.q.tau)-numel(p.tau),'complete',node.actual.passed); %#ok<AGROW>
   if node.actual.passed
    completeCount=completeCount+1;
    candidate=struct('q',node.q,'actual',node.actual,'found_s',toc(searchClock));
    candidates{end+1}=candidate; %#ok<AGROW>
    [~,ord]=sort(cellfun(@(n)n.actual.total_dv_km_s,candidates));
    candidates=candidates(ord(1:min(3,numel(ord))));
    if isequal(candidates{1}.q,node.q)
     fprintf('COMPLETE depth=%d t=%.1f J=%.9f (screened)\n',depth,toc(searchClock),node.actual.total_dv_km_s);
    end
    continue
   end
   % Safe only within this append-only diagnostic: future costs nonnegative.
   % Do NOT use a screened candidate to tighten this physical cost bound.
   if node.actual.total_dv_km_s>=baseline.total_dv_km_s
    boundPruned=boundPruned+1; continue
   end
   if node.zero_gain>c.max_zero_gain||numel(node.q.tau)-numel(p.tau)>=max(2*depth,depth+4), continue; end
   open{end+1}=node; %#ok<AGROW>
   scores(end+1)=node.actual.total_dv_km_s+node.actual.inclination_penalty ...
    +c.time_weight*35*max(0,node.q.T/864000-node.actual.visit_count/35); %#ok<AGROW>
  end
  % Lowest-score child is popped first; siblings remain for backtracking.
  [~,ord]=sort(scores,'descend'); stack=[stack,open(ord)]; %#ok<AGROW>
 end
 searchSeconds=toc(searchClock); verificationClock=tic; checked={}; winner=[];
 for j=1:numel(candidates)
  cand=candidates{j}; [v,~]=ctocscreen.v4.replay(cand.q,eph,c,true);
  cand.independent=v; checked{end+1}=cand; %#ok<AGROW>
  if v.passed&&(isempty(winner)||v.total_dv_km_s<winner.independent.total_dv_km_s), winner=cand; end
 end
 g=struct('depth',depth,'cut_burn',first,'cut_day',p.T/86400,'prefix_q',p, ...
  'retained_visits',pa.visit_count,'prefix_dv',pa.total_dv_km_s,'seed',927001, ...
  'remaining_targets',find(pa.distance_km>1).','expansions',expansions,'children',childrenCount, ...
  'duplicates',duplicates,'bound_pruned',boundPruned,'root_restarts',restarts, ...
  'guidance_failures',guidanceFailures,'complete_screened',completeCount, ...
  'best_partial_visits',bestPartial,'unexpanded_stack',numel(stack),'trials',trials, ...
  'checked',{checked},'best',winner,'search_seconds',searchSeconds, ...
  'verification_seconds',toc(verificationClock),'total_seconds',toc(groupClock));
 report.groups{end+1}=g;
 if isempty(winner)
  fprintf('RESULT depth=%d NO independent complete bestPartial=%d search=%.1fs\n',depth,bestPartial,searchSeconds);
 else
  v=winner.independent;
  fprintf('RESULT depth=%d independent=%d visits=%d J=%.12f saving=%.9f maxDist=%.6g height=%.6f T=%.6fd search=%.1fs\n', ...
   depth,v.passed,v.visit_count,v.total_dv_km_s,baseline.total_dv_km_s-v.total_dv_km_s, ...
   max(v.distance_km),v.min_altitude_lower_km,winner.q.T/86400,searchSeconds);
 end
 save(fullfile(outdir,'report.mat'),'report');
end
report.signature_unchanged=isequal(signature,ctocscreen.v4.signature());
report.source_unchanged=strcmp(report.source_sha256,ctocscreen.v3FileHash(source));
assert(report.signature_unchanged&&report.source_unchanged);
report.finished=true; report.total_seconds=toc(totalClock);
save(fullfile(outdir,'report.mat'),'report');
fprintf('DONE %.1fs %s\n',report.total_seconds,outdir);
end
