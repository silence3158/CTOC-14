function s=coneStep(p,z,e,trust,lambda,fuel,budget)
%CONESTEP Linearized nonlinear constraints; exact pulse and encounter cones.
n=p.n; ne=numel(e.eq); ng=numel(e.g); M=p.M; K=p.K;
is=n+(1:M); ip=n+M+(1:ne); im=n+M+ne+(1:ne);
iv=n+M+2*ne+(1:K); ig=n+M+2*ne+K+(1:ng); ie=n+M+2*ne+K+ng+1; ns=ie;
f=zeros(ns,1); f([ip im iv ig ie])=lambda; if ~fuel, f([ip im iv ig ie])=1; end
if fuel, f(is)=1; end
Ae=sparse(ne,ns); Ae(:,1:n)=e.Jeq; Ae(:,ip)=-speye(ne); Ae(:,im)=speye(ne);
Ai=sparse(ng+size(p.A,1),ns); Ai(1:ng,1:n)=e.Jg; Ai(1:ng,ig)=-speye(ng);
Ai(ng+1:end,1:n)=p.A;
bi=[-e.g;-p.A*z]; be=-e.eq;
lb=zeros(ns,1); ub=inf(ns,1);
lb(1:n)=max(-trust,p.lb-z); ub(1:n)=min(trust,p.ub-z);
soc=[];
for k=1:M
 A=sparse(3,ns); A(:,p.iu(:,k))=eye(3); d=zeros(ns,1); d(is(k))=1;
 item=secondordercone(A,-z(p.iu(:,k)),d,0);
 if isempty(soc), soc=item; else, soc(end+1)=item; end %#ok<AGROW>
end
for j=1:K
 A=sparse(3,ns); A(:,1:n)=e.Jvisit(3*(j-1)+(1:3),:); d=zeros(ns,1); d(iv(j))=1;
 item=secondordercone(A,-e.visit(:,j),d,-e.radius);
 if isempty(soc), soc=item; else, soc(end+1)=item; end %#ok<AGROW>
end
A=sparse(3,ns); A(:,1:n)=e.Jecc; d=zeros(ns,1); d(ie)=1;
item=secondordercone(A,-e.ecc,d,-e.ecc_radius);
if isempty(soc), soc=item; else, soc(end+1)=item; end
options=optimoptions('coneprog','Display','off','MaxIterations',p.config.solver_iterations, ...
 'MaxTime',max(.01,min(budget,p.config.solver_seconds)), ...
 'ConstraintTolerance',p.config.solver_tolerance,'OptimalityTolerance',p.config.solver_tolerance);
clock=tic; [w,value,flag,output]=coneprog(f,soc,Ai,bi,Ae,be,lb,ub,options);
s=struct('ok',false,'d',[],'exitflag',flag,'output',output,'seconds',toc(clock), ...
 'model_value',value,'slack',Inf,'nvar',ns,'neq',ne,'ncone',numel(soc));
if isempty(w)||any(~isfinite(w)), return; end
linear=max([0;Ai*w-bi;abs(Ae*w-be);lb-w;w-ub]); conic=0;
for j=1:numel(soc)
 conic=max(conic,norm(soc(j).A*w-soc(j).b)-(soc(j).d.'*w-soc(j).gamma));
end
s.feasibility=max(linear,conic); s.ok=s.feasibility<1e-5;
if s.ok
 s.d=w(1:n); s.slack=sum(w([ip im iv ig ie]));
end
end
