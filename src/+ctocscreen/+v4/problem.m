function p=problem(candidate,ids,theta,scope,focus,eph,c)
%PROBLEM Fixed event chart with affine auxiliary-node times and free controls.
% scope 'tail': focus is a time; every burn at or after it is open, the
% trajectory before the anchor t0 is left exactly as replayed (not integrated),
% and only encounters after t0 enter the subproblem. t0 is the later of the
% last fixed burn and the last active witness before the first open burn.
q=candidate.q; m=eph.model; ids=ids(:); theta=theta(:);
assert(numel(ids)==numel(theta)&&numel(unique(ids))==numel(ids));
assert(all(isfinite(theta))&&all(theta>=0&theta<=q.T));
levels={'local','expanded','full','tail'}; level=find(strcmp(scope,levels)); assert(~isempty(level));
tail=level==4; t0=0; bi=(1:numel(q.tau)).';
if tail
 k0=find(q.tau>=focus,1); assert(~isempty(k0),'ctocscreen:v4:tailEmpty','No burn in the tail window.');
 bi=(k0:numel(q.tau)).'; prior=0; if k0>1, prior=q.tau(k0-1); end
 t0=max([prior;theta(theta<q.tau(k0))]);
 keep=theta>t0; ids=ids(keep); theta=theta(keep);
end
tauP=q.tau(bi); uP=q.u(bi,:);
M=numel(bi); K=numel(ids); raw=[t0;tauP(:);theta;q.T]; rel=raw-t0;
[base,order]=sort(raw); nt=numel(raw)-1;
Ebase=sparse(numel(base),nt);
for j=1:numel(base), if order(j)>1, Ebase(j,order(j)-1)=1; end, end
Et=sparse(0,nt); baseNode=zeros(numel(base),1); times=[];
for j=1:numel(base)
 Et(end+1,:)=Ebase(j,:); times(end+1,1)=base(j); baseNode(j)=numel(times); %#ok<AGROW>
 if j<numel(base)
  segments=max(1,ceil((base(j+1)-base(j))/c.max_shooting_arc_s));
  for a=1:segments-1
   f=a/segments; Et(end+1,:)=(1-f)*Ebase(j,:)+f*Ebase(j+1,:); %#ok<AGROW>
   times(end+1,1)=(1-f)*base(j)+f*base(j+1); %#ok<AGROW>
  end
 end
end
N=numel(times); bn=zeros(M,1); wn=zeros(K,1);
for j=1:M, bn(j)=baseNode(order==1+j); end
for j=1:K, wn(j)=baseNode(order==1+M+j); end
L=m.re+600; V=sqrt(m.mu/L); st=L/V; sx=[L*ones(3,1);V*ones(3,1)];
ix=reshape(1:6*N,6,N); it=6*N+(1:nt); iu=reshape(6*N+nt+(1:3*M),3,M);
n=6*N+nt+3*M; zz=zeros(n,1);
X=ctocscreen.v4.stateAt(candidate.trace,times,'pre');
% The tail anchor is after every fixed burn, so it carries their effect.
if tail, X(:,1)=ctocscreen.v4.stateAt(candidate.trace,t0,'post'); end
% At coincident slots only preceding impulses have already occurred.
for j=1:N
 for k=find(tauP(:)==times(j)).'
  if bn(k)<j, X(4:6,j)=X(4:6,j)+uP(k,:).'; end
 end
end
zz(ix(:))=reshape(X./sx,[],1); zz(it)=rel(2:end)/st;
zz(iu(:))=reshape(uP.'/V,[],1);
lb=-inf(n,1); ub=inf(n,1); lb(it)=0; ub(it)=(m.horizon_s-t0)/st;
% Witness windows: local/expanded/full as before; tail uses the expanded radius.
radius=c.window_radius_s*2^(min(level,2)-1+(level==3));
for j=1:K
 col=it(M+j); lb(col)=max(0,(theta(j)-t0-radius)/st); ub(col)=min(m.horizon_s-t0,theta(j)-t0+radius)/st;
end
free=true(n,1); free(ix(:,1))=false; free(iu(:))=false; free(it(1:M))=false;
before=find(tauP<=focus); after=find(tauP>focus);
if tail
 open=1:M;
elseif level==3
 open=1:M; free(ix(:,1))=true;
else
 count=level; open=unique([before(max(1,end-count+1):end);after(1:min(count,numel(after)))]);
 if isempty(before)||M<=level, free(ix(:,1))=true; end
end
free(iu(:,open))=true; free(it(open))=true;
lb(~free)=zz(~free); ub(~free)=zz(~free);
A=sparse(N-1,n); A(:,it)=Et(1:end-1,:)-Et(2:end,:);
fixedCost=sum(vecnorm(q.u(setdiff(1:numel(q.tau),bi),:),2,2));
p=struct('q',q,'ids',ids,'scope',scope,'focus',focus,'N',N,'M',M,'K',K,'n',n, ...
 'ix',ix,'it',it,'iu',iu,'Et',Et,'bn',bn,'wn',wn,'sx',sx,'sv',V,'st',st, ...
 'res_scale',sx,'length_scale',L,'z0',zz,'lb',lb,'ub',ub,'A',A, ...
 'height_fractions',[.25 .5 .75],'free',free,'opened_burns',bi(open).','model',m,'eph',eph,'config',c, ...
 'raw_times',rel(2:end),'t0',t0,'tail',tail,'burn_index',bi,'fixed_cost',fixedCost);
end
