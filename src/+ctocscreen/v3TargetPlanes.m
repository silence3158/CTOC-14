function [planes,counts,members]=v3TargetPlanes(eph,c)
%V3TARGETPLANES Distinct target orbital planes with a representative target each.
% Search guidance for root generation: the initial plane is a free design
% variable, so the deterministic root families should cover the planes that
% actually appear in the target set. Not a physical constraint, not a result.
if nargin<2, c=ctocscreen.v3Defaults(); end
x=eph.states0; h=cross(x(:,1:3),x(:,4:6),2);
inc=acosd(h(:,3)./vecnorm(h,2,2));
edges=0:c.target_plane_bin_deg:180;
bin=discretize(inc,edges);
planes=cell(0,1); counts=[]; members={};
for b=unique(bin(~isnan(bin))).'
 ids=find(bin==b);
 % planes{k} is the median true inclination of the bin, so the representative
 % plane sits on the population rather than on a bin edge; members{k} lists the
 % targets belonging to that plane and is used to pick a representative target.
 planes{end+1,1}=median(inc(ids)); %#ok<AGROW>
 counts(end+1,1)=numel(ids); %#ok<AGROW>
 members{end+1}=ids(:).'; %#ok<AGROW>
end
[counts,order]=sort(counts,'descend');
planes=planes(order); members=members(order);
end

