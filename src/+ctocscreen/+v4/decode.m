function [q,X,t]=decode(p,z)
%DECODE No projection, re-aiming, or state reset.
% Frozen controls are copied from the task, not rescaled: (u/V)*V is not
% bit-exact, and exact frozen values keep physical keys and prefix reuse valid.
% Times are clamped monotone: coneprog meets the linear ordering only to tolerance.
% Tail tasks shift times by t0 and replace only the burns in p.burn_index.
q=p.q; X=reshape(z(p.ix(:)),6,p.N).*p.sx;
values=p.st*z(p.it); fixed=~p.free(p.it); values(fixed)=p.raw_times(fixed);
t=p.t0+min(p.model.horizon_s-p.t0,cummax(max(0,p.Et*values)));
if ~p.tail&&p.free(p.ix(1,1)), q.x0=X(:,1); end
u=reshape(z(p.iu(:)),3,p.M).'*p.sv; frozen=~p.free(p.iu(1,:));
u(frozen,:)=p.q.u(p.burn_index(frozen),:);
q.tau(p.burn_index)=t(p.bn); q.u(p.burn_index,:)=u; q.T=t(end);
q.witness(p.ids)=t(p.wn);
end
