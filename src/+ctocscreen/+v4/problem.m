function p=problem(candidate,ids,theta,scope,focus,eph,c)
%PROBLEM Fixed event chart with affine auxiliary-node times and free controls.
q=candidate.q; m=eph.model; ids=ids(:); theta=theta(:);
assert(numel(ids)==numel(theta)&&numel(unique(ids))==numel(ids));
assert(all(isfinite(theta))&&all(theta>=0&theta<=q.T));
M=numel(q.tau); K=numel(ids); raw=[0;q.tau(:);theta;q.T];
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
% At coincident slots only preceding impulses have already occurred.
for j=1:N
 for k=find(q.tau(:)==times(j)).'
  if bn(k)<j, X(4:6,j)=X(4:6,j)+q.u(k,:).'; end
 end
end
zz(ix(:))=reshape(X./sx,[],1); zz(it)=raw(2:end)/st;
zz(iu(:))=reshape(q.u.'/V,[],1);
lb=-inf(n,1); ub=inf(n,1); lb(it)=0; ub(it)=m.horizon_s/st;
levels={'local','expanded','full'}; level=find(strcmp(scope,levels)); assert(~isempty(level));
radius=c.window_radius_s*2^(level-1);
for j=1:K
 col=it(M+j); lb(col)=max(0,(theta(j)-radius)/st); ub(col)=min(m.horizon_s,(theta(j)+radius))/st;
end
free=true(n,1); free(ix(:,1))=false; free(iu(:))=false; free(it(1:M))=false;
before=find(q.tau<=focus); after=find(q.tau>focus);
if level==3
 open=1:M; free(ix(:,1))=true;
else
 count=level; open=unique([before(max(1,end-count+1):end);after(1:min(count,numel(after)))]);
 if isempty(before)||M<=level, free(ix(:,1))=true; end
end
free(iu(:,open))=true; free(it(open))=true;
lb(~free)=zz(~free); ub(~free)=zz(~free);
A=sparse(N-1,n); A(:,it)=Et(1:end-1,:)-Et(2:end,:);
p=struct('q',q,'ids',ids,'scope',scope,'focus',focus,'N',N,'M',M,'K',K,'n',n, ...
 'ix',ix,'it',it,'iu',iu,'Et',Et,'bn',bn,'wn',wn,'sx',sx,'sv',V,'st',st, ...
 'res_scale',sx,'length_scale',L,'z0',zz,'lb',lb,'ub',ub,'A',A, ...
 'height_fractions',[.25 .5 .75],'free',free,'opened_burns',open,'model',m,'eph',eph,'config',c, ...
 'raw_times',raw(2:end));
end
