function [beam,estimates]=selectBeam(nodes,c,eph,stream)
%SELECTBEAM Layered beam: compare candidates only within equal visit counts.
% Rank within a layer by f = J + H + P (km/s): J raw spent Delta-V, H the
% Lambert crossing estimate of remaining cost (estimate.m, heuristic), P a
% time-overrun penalty. Layers after Beam P-ACO depth levels. Within a layer
% the plane normal and orbit shape keep representatives apart.
if nargin<3, eph=[]; end
if nargin<4, stream=RandStream('mt19937ar','Seed',0); end
estimates=0;
valid={}; keys={};
for k=1:numel(nodes)
 n=nodes{k}; a=n.actual;
 if ~a.initial_passed||~a.height_passed||strcmp(a.status,'propagation_failure')||~isfinite(a.total_dv_km_s), continue; end
 if n.zero_gain>c.max_zero_gain, continue; end
 key=ctocscreen.v4.controlKey(n.q); match=find(strcmp(keys,key),1);
 if isempty(match)
  valid{end+1}=n; keys{end+1}=key; %#ok<AGROW>
 elseif n.attempts>valid{match}.attempts
  valid{match}=n;
 end
end
beam={}; if isempty(valid), return; end
for k=1:numel(valid)
 if ~isfield(valid{k},'heuristic_H')||isnan(valid{k}.heuristic_H)
  valid{k}.heuristic_H=NaN;
  if ~isempty(eph)
   if strcmp(c.heuristic_kind,'lambert')
    valid{k}.heuristic_H=ctocscreen.v4.estimate(valid{k},eph,c,stream);
   else
    [~,info]=ctocscreen.v4.events(valid{k},eph,c); valid{k}.heuristic_H=info.H;
   end
   estimates=estimates+1;
  end
 end
end
count=cellfun(@(n)n.actual.visit_count,valid); f=cellfun(@(n)merit(n,c),valid);
layers=unique(count);
% Deepest layers first; every nonempty layer keeps at least one representative.
quota=zeros(size(layers)); total=min(c.beam_width,numel(valid));
order=numel(layers):-1:1;
for k=order, if sum(quota)<total, quota(k)=1; end, end
while sum(quota)<total
 grown=false;
 for k=order
  if sum(quota)>=total, break; end
  if quota(k)<sum(count==layers(k)), quota(k)=quota(k)+1; grown=true; end
 end
 if ~grown, break; end
end
F=zeros(numel(valid),9);
for k=1:numel(valid)
 x=valid{k}.actual.final_state; h=cross(x(1:3),x(4:6)); h=h/max(norm(h),eps);
 R=norm(x(1:3)); en=dot(x(4:6),x(4:6))/2-398600.4415/R;
 F(k,:)=[h.',R/10000,en/5,valid{k}.q.T/86400,x(4:6).'/2];
end
selected=[];
for k=order
 members=find(count==layers(k)); if quota(k)==0, continue; end
 [~,rank]=sort(f(members)); pick=members(rank(1));
 while numel(pick)<quota(k)
  rest=setdiff(members,pick,'stable');
  if isempty(rest), break; end
  d=zeros(numel(rest),1);
  for j=1:numel(rest), d(j)=min(vecnorm(F(pick,:)-F(rest(j),:),2,2)); end
  % Prefer low merit, but reward distinct end geometry within the layer.
  score=reshape(f(rest),[],1)-c.layer_diversity*min(1,d);
  [~,j]=min(score); pick(end+1)=rest(j); %#ok<AGROW>
 end
 selected=[selected,pick]; %#ok<AGROW>
end
beam=valid(selected);
end
function f=merit(n,c)
J=n.actual.total_dv_km_s; H=n.heuristic_H; if isnan(H), H=0; end
k=n.actual.visit_count; used=n.q.T/864000;
over=max(0,used-(k/35+c.time_slack));
f=J+H+c.time_weight*over*35;
end
