function c=mutateCandidate(c,p,cfg,stream)
%MUTATECANDIDATE Mix order, time and all initial-orbit degrees of freedom.
if isfield(c,'branch_refs'), c=rmfield(c,'branch_refs'); end
c.branch_ids=zeros(1,35); c.parent_id=c.candidate_id;
switch randi(stream,4)
 case 1
  ij=randperm(stream,35,2); c.order(ij)=c.order(fliplr(ij));
 case 2
  ij=randperm(stream,35,2); item=c.order(ij(1)); c.order(ij(1))=[];
  c.order=[c.order(1:ij(2)-1) item c.order(ij(2):end)];
 case 3
  ix=randperm(stream,35,randi(stream,5));
  c.tof_s(ix)=max(cfg.min_tof_s,c.tof_s(ix).*exp(.15*randn(stream,size(ix))));
 otherwise
  q=c.initial_q+[2 .0001 .0001 .04 .08 .15].*randn(stream,1,6);
  q(1)=max(p.re_km+590,min(p.re_km+610,q(1)));
  e=hypot(q(2),q(3)); if e>=.001, q(2:3)=q(2:3)*(.000999/e); end
  q(4)=max(0,min(pi,q(4))); q(5:6)=mod(q(5:6),2*pi); c.initial_q=q;
  c.wait_s=max(0,c.wait_s+300*randn(stream));
end
if c.wait_s+sum(c.tof_s)>p.horizon_s
 budget=p.horizon_s-c.wait_s-35*cfg.min_tof_s;
 c.tof_s=cfg.min_tof_s+(c.tof_s-cfg.min_tof_s)*max(0,budget)/sum(c.tof_s-cfg.min_tof_s);
end
c.generation_method="mutation";
end
