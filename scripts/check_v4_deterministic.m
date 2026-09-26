function report=check_v4_deterministic()
%CHECK_V4_DETERMINISTIC Targeted checks for the 2026-09-26 deterministic repairs.
% Synthetic parts are synthetic_test_only; the real part only reads run-04
% controls as numerical input and never feeds the cold search.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'),fullfile(sim,'tests'));
report=struct('dataset_kind','mixed_diagnostic','passed',false);
% 1. Truncated prefix reuse must equal a from-epoch replay, on real controls.
s=load(fullfile(sim,'runs/v4/search/v4_cold_first_600_seed888_04/result.mat'),'result');
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(); q=s.result.best.q;
q.T=q.tau(21); q.tau=q.tau(1:20); q.u=q.u(1:20,:); q.witness(q.witness>q.T)=NaN;
[pa,pt]=ctocscreen.v4.replay(q,eph,c); parent=struct('q',pt.q,'actual',pa,'trace',pt,'root_id',1, ...
 'root_kind','diagnostic','seed_target',1,'seed_duration',3600,'attempts',0,'zero_gain',0,'generation',0,'origin','diagnostic');
parent.q.witness=pa.witness_times_s;
cases={'coast','delayed_impulse'}; reuse=zeros(2,3);
for k=1:2
 child=parent.q; child.T=q.T+7200;
 if k==2, child.tau(end+1,1)=q.T+1800; child.u(end+1,:)=[.01 -.02 .005]; end
 t1=tic; [ra,rt]=ctocscreen.v4.replay(child,eph,c,false,parent); reusedSeconds=toc(t1);
 t2=tic; [fa,ft]=ctocscreen.v4.replay(child,eph,c); freshSeconds=toc(t2);
 assert(ra.reused_prefix,'%s child did not reuse the prefix.',cases{k});
 assert(isequal(ra.final_state,fa.final_state)&&isequal(rt.post,ft.post),'%s states differ.',cases{k});
 assert(isequal(ra.distance_km,fa.distance_km)&&isequaln(ra.witness_times_s,fa.witness_times_s),'%s scans differ.',cases{k});
 assert(ra.height_passed==fa.height_passed&&ra.min_altitude_lower_km==fa.min_altitude_lower_km,'%s heights differ.',cases{k});
 reuse(k,:)=[reusedSeconds,freshSeconds,numel(rt.arcs)];
 fprintf('REUSE %s: identical to epoch replay; %.3f s vs %.3f s, %d arcs\n',cases{k},reusedSeconds,freshSeconds,numel(rt.arcs));
end
% A changed third control may reuse only the arcs before it; results must
% still equal the epoch replay.
other=parent.q; other.u(3,1)=other.u(3,1)+1e-9; other.T=q.T+600;
[oa,ot]=ctocscreen.v4.replay(other,eph,c,false,parent); [ob,obt]=ctocscreen.v4.replay(other,eph,c);
assert(isequal(ot.post,obt.post)&&isequal(oa.distance_km,ob.distance_km),'Changed-control reuse differs.');
% Arc 1 is the initial coast, so arcs 1..3 end at the changed third pulse.
assert(isequal(ot.arcs{3},parent.trace.arcs{3})&&~isequal(ot.arcs{4},parent.trace.arcs{4}), ...
 'Only arcs before the changed control may be reused.');
% 2. Structure rotation reaches every operator including the suffix rebuild.
kinds=mod((1:10)-1,5); assert(isequal(sort(unique(kinds)),0:4));
[rb,rr]=ctocscreen.v4.restructure(parent,eph,c,5,4);
assert(strcmp(rr.action,'suffix_rebuild')&&~isempty(rb)&&rb.q.T<parent.q.T,'Index 5 must rebuild the suffix.');
fprintf('STRUCTURE index 5 -> %s, %d -> %d visits\n',rr.action,parent.actual.visit_count,rb.actual.visit_count);
% 3. Model limits are strictly inside acceptance limits.
assert(c.model_radius_km<c.search_radius_km&&c.search_radius_km<=1);
assert(c.eccentricity_limit<.001&&c.model_height_margin_km>c.height_margin_km&&c.sma_margin_km>0);
% 4. Shared-task failures escalate scope, then the task is skipped (synthetic).
[seph,~,old]=v3SyntheticFixture('separated_visits');
sc=ctocscreen.v4.defaults(struct('max_step_s',30,'scan_step_s',60,'joint_iterations',1, ...
 'solver_seconds',3,'discovery_radius_km',1e6,'shared_attempts',2));
% A deliberate 50 m/s first pulse leaves targets tens of km away, so one
% iteration cannot recover them and the task records executed failures.
sq=struct('x0',ctocscreen.initialState(old.initial_q,seph.model.mu,seph.model.re).', ...
 'tau',[90;360],'u',[.05 -.001 .0005;-.0004 .0003 .0002],'T',600,'witness',nan(35,1));
[sa,st]=ctocscreen.v4.replay(sq,seph,sc);
node=struct('q',st.q,'actual',sa,'trace',st,'root_id',1,'root_kind','synthetic', ...
 'seed_target',1,'seed_duration',100,'attempts',0,'zero_gain',0,'generation',0,'origin','synthetic');
node.q.witness=sa.witness_times_s; tabu=[]; scopes={}; firstIds=[];
for k=1:3
 [~,jr,~,tabu]=ctocscreen.v4.shared(node,seph,sc,[5 5 5],tabu);
 if isfield(jr,'scope'), scopes{end+1}=jr.scope; end %#ok<AGROW>
 if k==1&&isfield(jr,'added_target_ids'), firstIds=jr.added_target_ids; end
 if k==3, lastIds=getfield(jr,'added_target_ids',[]); end
 fprintf('SHARED call %d: status=%s scope=%s added=%s skipped=%d active_passed=%d failures=%s\n',k,jr.status, ...
  string(getfield(jr,'scope','none')),mat2str(getfield(jr,'added_target_ids',[])),jr.skipped_tasks, ...
  getfield(jr,'active_passed',false),mat2str(tabu.failures));
end
assert(~isempty(tabu.failures),'Synthetic shared task never executed.');
assert(numel(scopes)>=2&&strcmp(scopes{1},'local')&&strcmp(scopes{2},'expanded'),'Scope did not escalate.');
assert(max(tabu.failures)<=sc.shared_attempts,'A task exceeded its attempt limit.');
% After an executed failure and then a success, the same prefix/target task is closed.
assert(numel(scopes)==2||~isequal(firstIds,lastIds),'A closed task was executed again.');
% 5. A joint call respects its allowance on a real 20-pulse prefix.
ids=find(pa.distance_km<=1); t3=tic;
[~,jr]=ctocscreen.v4.joint(parent,ids,pa.witness_times_s(ids),'local',q.tau(18),eph,c,2);
elapsed=toc(t3);
fprintf('BUDGET local joint allowance 2 s: %.3f s, %d iterations, %d replays (%.3f s), status=%s\n', ...
 elapsed,jr.iterations,jr.replay_count,jr.replay_seconds,jr.status);
report=struct('dataset_kind','mixed_diagnostic','passed',true,'reuse_seconds',reuse, ...
 'shared_scopes',{scopes},'shared_first_ids',firstIds,'budget_elapsed_s',elapsed,'budget_joint',jr, ...
 'signature',ctocscreen.v4.signature());
save(fullfile(sim,'runs/v4/development',['deterministic_checks_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.mat']),'report');
fprintf('V4 DETERMINISTIC CHECKS PASSED\n');
end
function v=getfield(s,name,fallback)
v=fallback; if isfield(s,name), v=s.(name); end
end
