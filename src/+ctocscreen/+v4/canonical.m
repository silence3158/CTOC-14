function q=canonical(q,m)
%CANONICAL Explicit control normalization outside the smooth evaluator.
q.x0=q.x0(:); q.tau=q.tau(:);
assert(numel(q.x0)==6&&all(isfinite(q.x0)));
assert(size(q.u,1)==numel(q.tau)&&size(q.u,2)==3);
assert(all(isfinite(q.u(:)))&&all(isfinite(q.tau)));
assert(isscalar(q.T)&&isfinite(q.T)&&q.T>=0&&q.T<=m.horizon_s);
assert(all(q.tau>=0&q.tau<=q.T),'ctocscreen:v4:time','Pulse outside mission.');
[q.tau,order]=sort(q.tau); q.u=q.u(order,:);
if ~isempty(q.tau)
 [tt,~,group]=unique(q.tau); uu=zeros(numel(tt),3);
 for k=1:numel(group), uu(group(k),:)=uu(group(k),:)+q.u(k,:); end
 q.tau=tt; q.u=uu;
 % Exact zero pulses have no physical effect. Near-zero pulses are retained.
 keep=any(q.u~=0,2)&q.tau<q.T;
 q.tau=q.tau(keep); q.u=q.u(keep,:);
end
if ~isfield(q,'witness'), q.witness=nan(35,1); end
q.witness=q.witness(:); assert(numel(q.witness)==35);
q.schema_version='free_trajectory_v4_1'; q.dynamics_id='central_j2';
end
