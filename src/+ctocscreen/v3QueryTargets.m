function [r,v,a] = v3QueryTargets(cache,ids,t,mode)
%V3QUERYTARGETS Fast target states from a loaded cache (no disk access).
% 'pairs' (default): scalar expansion or equal lengths; outputs N-by-3.
% 'grid': all IDs at all times; outputs numel(ids)-by-3-by-numel(t).
% [] IDs means 1:35. Velocity is the exact derivative of position polynomial.
% Acceleration is evaluated from the J2 force model, not polynomial curvature.
if nargin<4, mode='pairs'; end
if isempty(ids), ids=1:35; end
validateattributes(ids,{'double'},{'vector','integer','>=',1,'<=',35,'finite','real'});
validateattributes(t,{'double'},{'vector','nonempty','real','finite'});
assert(strcmp(cache.schema_version,'target_ephemeris_v3_1'),'ctocscreen:v3:schema','Unsupported cache.');
assert(all(t>=0 & t<=cache.model.horizon_s),'ctocscreen:v3:timeRange','No extrapolation allowed.');
ids=ids(:); t=t(:); ni=numel(ids); nt=numel(t);
if strcmp(mode,'grid')
 ids=repmat(ids,nt,1); t=repelem(t,ni);
elseif strcmp(mode,'pairs')
 if isscalar(ids), ids=repmat(ids,nt,1); end
 if isscalar(t), t=repmat(t,numel(ids),1); end
 assert(numel(ids)==numel(t),'ctocscreen:v3:shape','Use grid or equal-length pairs.');
else
 error('ctocscreen:v3:mode','Mode must be pairs or grid.');
end
r=zeros(numel(t),3); v=r;
for id=unique(ids).'
 sel=find(ids==id); block=cache.targets{id}; h=block.step_s;
 idx=min(floor(t(sel)/h)+1,size(block.coef,1));
 u=(t(sel)-(idx-1)*h)/h; c=block.coef(idx,:,:);
 rp=c(:,:,6); vp=5*c(:,:,6);
 for d=5:-1:1
  rp=rp.*u+c(:,:,d);
  if d>1, vp=vp.*u+(d-1)*c(:,:,d); end
 end
 r(sel,:)=rp; v(sel,:)=vp/h;
end
if nargout>2, a=ctocscreen.v3J2Acceleration(t,r,cache.model); end
if strcmp(mode,'grid')
 r=permute(reshape(r,ni,nt,3),[1 3 2]);
 v=permute(reshape(v,ni,nt,3),[1 3 2]);
 if nargout>2, a=permute(reshape(a,ni,nt,3),[1 3 2]); end
end
end
