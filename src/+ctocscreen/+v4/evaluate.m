function e=evaluate(p,z,derivatives)
%EVALUATE Exact nonlinear event constraints and sparse analytic Jacobian.
if nargin<3, derivatives=true; end
[q,X,t]=ctocscreen.v4.decode(p,z); c=p.config; m=p.model;
assert(all(t>=0)&all(diff(t)>=0)&t(end)<=m.horizon_s,'ctocscreen:v4:eventOrder','Invalid event chart.');
n=p.n; N=p.N; M=p.M; K=p.K; nf=numel(p.height_fractions);
e.eq=zeros(6*(N-1),1); e.Jeq=sparse(numel(e.eq),n);
e.g=zeros(2+N+nf*(N-1),1); e.Jg=sparse(numel(e.g),n);
e.visit=zeros(3,K); e.Jvisit=sparse(3*K,n);
e.radius=c.search_radius_km; e.ecc_radius=c.eccentricity_limit/.001;
o=ctocscreen.v4.initialOrbit(q.x0,m); width=1/(m.re+590)-1/(m.re+610);
e.g(1:2)=[1/(m.re+610)-o.alpha;o.alpha-1/(m.re+590)]/width;
e.Jg(1:2,p.ix(:,1))=[-o.alpha_jac;o.alpha_jac].*p.sx.'/width;
e.ecc=o.evec/.001; e.Jecc=sparse(3,n);
e.Jecc(:,p.ix(:,1))=o.ecc_jac.*p.sx.'/.001;
U=zeros(3,N); for k=1:M, U(:,p.bn(k))=q.u(k,:).'; end
sx=p.sx; rs=p.res_scale;
for l=1:N-1
 x=X(:,l)+[zeros(3,1);U(:,l)];
 [yf,P,sol]=ctocscreen.v3Arc(x,t(l),t(l+1),m,c,derivatives);
 rows=6*(l-1)+(1:6); e.eq(rows)=(X(:,l+1)-yf)./rs;
 fa=[x(4:6);ctocscreen.v3Force(t(l),x(1:3),m)];
 fb=[yf(4:6);ctocscreen.v3Force(t(l+1),yf(1:3),m)];
 b=find(p.bn==l);
 if derivatives
  e.Jeq(rows,p.ix(:,l))=-P.*sx.'./rs;
  e.Jeq(rows,p.ix(:,l+1))=diag(sx./rs);
  if ~isempty(b), e.Jeq(rows,p.iu(:,b))=-P(:,4:6)*p.sv./rs; end
  e.Jeq(rows,p.it)=(P*fa./rs)*p.st*p.Et(l,:)-(fb./rs)*p.st*p.Et(l+1,:);
 end
 for f=1:nf
  a=p.height_fractions(f); tm=(1-a)*t(l)+a*t(l+1);
  if isempty(sol), ym=x; Pm=eye(6);
  else
   yy=deval(sol,tm); ym=yy(1:6);
   if derivatives, Pm=reshape(yy(7:end),6,6); end
  end
  row=2+N+(l-1)*nf+f; rr=norm(ym(1:3)); e.g(row)=m.re+200+c.height_margin_km-rr;
  if derivatives
   h=[-ym(1:3).'/rr,zeros(1,3)]; fm=[ym(4:6);ctocscreen.v3Force(tm,ym(1:3),m)];
   e.Jg(row,p.ix(:,l))=h*Pm.*sx.';
   if ~isempty(b), e.Jg(row,p.iu(:,b))=h*Pm(:,4:6)*p.sv; end
   e.Jg(row,p.it)=h*((1-a)*fm-Pm*fa)*p.st*p.Et(l,:)+h*(a*fm)*p.st*p.Et(l+1,:);
  end
 end
end
for l=1:N
 rr=norm(X(1:3,l)); e.g(2+l)=m.re+200+c.height_margin_km-rr;
 e.Jg(2+l,p.ix(1:3,l))=-X(1:3,l).'/rr.*sx(1:3).';
end
for j=1:K
 l=p.wn(j); [rt,vt]=ctocscreen.v3QueryTargets(p.eph,p.ids(j),t(l));
 e.visit(:,j)=X(1:3,l)-rt.'; rows=3*(j-1)+(1:3);
 e.Jvisit(rows,p.ix(1:3,l))=diag(sx(1:3));
 e.Jvisit(rows,p.it)=-vt.'*p.st*p.Et(l,:);
end
e.J=sum(vecnorm(q.u,2,2)); e.V=violation(e);
e.max_position_defect=0; e.max_velocity_defect=0;
if N>1
 dd=reshape(e.eq,6,[]).*rs;
 e.max_position_defect=max(vecnorm(dd(1:3,:),2,1));
 e.max_velocity_defect=max(vecnorm(dd(4:6,:),2,1));
end
e.near_feasible=e.V<c.restore_tolerance&&e.max_position_defect<c.defect_position_km ...
 &&e.max_velocity_defect<c.defect_velocity_km_s;
assert(all(isfinite([e.eq;e.g;e.ecc;e.visit(:);e.J;e.V])),'ctocscreen:v4:nonfinite','Nonfinite constraints.');
end
function v=violation(e)
v=sum(abs(e.eq))+sum(max(0,e.g))+sum(max(0,vecnorm(e.visit,2,1)-e.radius)) ...
 +max(0,norm(e.ecc)-e.ecc_radius);
end
