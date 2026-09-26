function report=eval_v4_polish(file,budget)
%EVAL_V4_POLISH Does B still pay on an already complete trajectory?
% Runs full-history B (joint.m) repeatedly on one independently complete
% construction result, restarting from each improved incumbent, within a
% budget. Reports actual independent Delta-V before and after.
if nargin<2, budget=120; end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(); s=load(fullfile(sim,file));
if isfield(s,'result'), q=s.result.best.q; elseif isfield(s.report,'best'), q=s.report.best.q; else, q=s.report.q; end
[a,tr]=ctocscreen.v4.replay(q,eph,c);
node=struct('q',tr.q,'actual',a,'trace',tr,'root_id',1,'root_kind','evaluation','seed_target',1, ...
 'seed_duration',3600,'attempts',0,'zero_gain',0,'generation',0,'origin','evaluation','heuristic_H',NaN);
node.q.witness=a.witness_times_s; J0=a.total_dv_km_s; clock=tic; calls=0; trace=J0; warm=[];
while toc(clock)<budget
 ids=find(node.actual.distance_km<=1); calls=calls+1;
 [child,jr,warm]=ctocscreen.v4.joint(node,ids,node.actual.witness_times_s(ids),'full',node.q.T,eph,c, ...
  min(c.scope_seconds(3)*2,budget-toc(clock)),resumeIf(warm,node));
 if ~isempty(child)&&jr.active_passed&&child.actual.total_dv_km_s<node.actual.total_dv_km_s-1e-7
  child.q.witness=child.actual.witness_times_s; node=child; warm=[];
 end
 trace(end+1)=node.actual.total_dv_km_s; %#ok<AGROW>
 fprintf('POLISH call %d t=%.1f status=%s it=%d J=%.6f\n',calls,toc(clock),jr.status,jr.iterations,node.actual.total_dv_km_s);
end
[v,~]=ctocscreen.v4.replay(node.q,eph,c,true);
fprintf('POLISH %s: %.6f -> %.6f km/s (independent %d/35 passed=%d), %d calls, %.1f s\n', ...
 file,J0,v.total_dv_km_s,v.visit_count,v.passed,calls,toc(clock));
report=struct('evaluation_only',true,'file',file,'J0',J0,'J1',v.total_dv_km_s,'passed',v.passed, ...
 'trace',trace,'calls',calls);
end
function r=resumeIf(warm,node)
r=[];
if ~isempty(warm)&&isfield(warm,'physical_key')&&strcmp(warm.physical_key,ctocscreen.v4.controlKey(node.q))&&warm.resumable
 r=warm;
end
end
