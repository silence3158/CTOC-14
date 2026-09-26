function x=stateAt(tr,t,side)
%STATEAT State from a fixed-control trace, with an explicit pulse side.
if nargin<3, side='pre'; end
t=t(:).'; x=zeros(6,numel(t));
for k=1:numel(t)
 j=find(tr.times==t(k),1);
 if ~isempty(j)
  if strcmp(side,'pre'), x(:,k)=tr.pre(:,j); else, x(:,k)=tr.post(:,j); end
 else
  j=find(tr.arc_start<t(k)&tr.arc_end>t(k),1);
  assert(~isempty(j),'ctocscreen:v4:traceRange','Time outside propagated trace.');
  yy=deval(tr.arcs{j},t(k)); x(:,k)=yy(1:6);
 end
end
end
