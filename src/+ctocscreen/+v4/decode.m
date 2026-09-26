function [q,X,t]=decode(p,z)
%DECODE No projection, re-aiming, or state reset.
q=p.q; X=reshape(z(p.ix(:)),6,p.N).*p.sx;
t=p.st*(p.Et*z(p.it)); q.x0=X(:,1);
q.tau=t(p.bn); q.u=reshape(z(p.iu(:)),3,p.M).'*p.sv; q.T=t(end);
q.witness(p.ids)=t(p.wn);
end
