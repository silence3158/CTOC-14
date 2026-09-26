function [child,report]=rebuild(parent,eph,c,attempt)
%REBUILD Revoke suffix witnesses and resume from an earlier real control state.
clock=tic; child=[]; q=parent.q; M=numel(q.tau);
report=struct('status','no_rebuild','action','suffix_rebuild','iterations',0,'seconds',0, ...
 'parent_visits',parent.actual.visit_count,'retained_visits',0,'cut_time_s',NaN,'revoked_target_ids',[]);
if M==0, return; end
cost=vecnorm(q.u,2,2); [~,expensive]=max(cost);
over=find(cumsum(cost)>c.search_max_dv_km_s,1); if isempty(over), over=expensive; end
choices=unique([max(1,M-3),expensive,over,max(1,ceil(M/2))],'stable');
k=choices(1+mod(attempt-1,numel(choices))); cut=q.tau(k);
q.T=cut; q.tau=q.tau(1:k-1); q.u=q.u(1:k-1,:); q.witness(q.witness>cut)=NaN;
[a,tr]=ctocscreen.v4.replay(q,eph,c);
if ~a.initial_passed||~a.height_passed||strcmp(a.status,'propagation_failure')
 report.status='rebuild_replay_failed'; report.seconds=toc(clock); return
end
child=parent; child.q=tr.q; child.q.witness=a.witness_times_s; child.actual=a; child.trace=tr;
child.origin='suffix_rebuild'; child.root_kind='reconstructed_current_run';
child.attempts=0; child.zero_gain=0; child.generation=parent.generation+1; child.heuristic_H=NaN;
report.status='actual_prefix'; report.cut_time_s=cut; report.retained_visits=a.visit_count;
report.revoked_target_ids=find(parent.actual.distance_km<=1&a.distance_km>1).';
report.seconds=toc(clock);
end
