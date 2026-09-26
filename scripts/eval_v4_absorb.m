function report=eval_v4_absorb(file,label)
%EVAL_V4_ABSORB Step-3 check: how often can a following encounter be folded in?
% Walks the prefixes of a saved chain (diagnostic input only). At each prefix
% ending in an arrival, runs absorb(); records candidates, prescreen outcome,
% tail-B result, and compares the realized extra cost with the cheapest
% standalone Lambert leg to the same target from that node.
if nargin<2, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(); stream=RandStream('mt19937ar','Seed',5);
s=load(fullfile(sim,file)); Q=s.report.q; rows=zeros(0,9); clock=tic;
for k=2:numel(Q.tau)
 q=Q; q.T=Q.tau(k); q.tau=Q.tau(1:k-1); q.u=Q.u(1:k-1,:); q.witness(q.witness>q.T)=NaN;
 [a,tr]=ctocscreen.v4.replay(q,eph,c);
 node=struct('q',tr.q,'actual',a,'trace',tr,'root_id',1,'root_kind','e','seed_target',1,'seed_duration',3600, ...
  'attempts',0,'zero_gain',0,'generation',0,'origin','e','heuristic_H',NaN); node.q.witness=a.witness_times_s;
 [child,r]=ctocscreen.v4.absorb(node,eph,c,c.absorb_seconds);
 standalone=NaN; extra=NaN;
 if ~isnan(r.target)
  [~,d]=ctocscreen.v4.estimate(node,eph,c,stream); standalone=d.per_target(r.target);
 end
 if ~isempty(child), extra=child.actual.total_dv_km_s-a.total_dv_km_s; end
 rows(end+1,:)=[k-1,a.visit_count,r.candidates,r.screened,r.predicted_extra,r.prescreen_residual_km,extra,standalone,r.seconds]; %#ok<AGROW>
 fprintf('ABSORB prefix %2d visits %2d: cand %d status %-22s target %2d pred %+.3f actual %+.3f standalone %.3f (%.1f s)\n', ...
  k-1,a.visit_count,r.candidates,r.status,r.target,r.predicted_extra,extra,standalone,r.seconds);
end
done=isfinite(rows(:,7));
fprintf('ABSORB summary: prefixes %d, with candidates %d, passed prescreen %d, absorbed %d; median actual extra %.3f vs standalone %.3f; %.1f s\n', ...
 size(rows,1),sum(rows(:,3)>0),sum(isfinite(rows(:,6))),sum(done),median(rows(done,7)),median(rows(done,8)),toc(clock));
report=struct('diagnostic_only',true,'file',file,'rows',rows,'signature',ctocscreen.v4.signature());
save(fullfile(sim,'runs/v4/development',['absorb_eval_' label '.mat']),'report');
end
