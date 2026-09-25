function beam=v3SelectBeam(nodes,width,c)
%V3SELECTBEAM Merge near states before coverage/coast quotas and diversity fill.
if nargin<3, c=ctocscreen.v3Defaults(); end
if isempty(nodes), beam={}; return; end
score=cellfun(@(n)n.estimate,nodes); [~,order]=sort(score);
kept=[];
for k=order
 if isempty(kept)||all(cellfun(@(n)distance(nodes{k},n,c)>=1,nodes(kept)))
  kept(end+1)=k; %#ok<AGROW>
 end
end
nodes=nodes(kept); score=score(kept); count=min(width,numel(nodes));
coverage=cellfun(@(n)sum(n.visited),nodes);
over=false(size(coverage));
if c.cost_guidance_enabled, over=cellfun(@(n)n.J>c.search_max_dv_km_s,nodes); end
[~,progress]=sortrows([over(:),-coverage(:),score(:)]);
chosen=progress(1); progressQuota=min(count,max(1,ceil(width/2)));
while numel(chosen)<progressQuota
 left=setdiff(progress(:).',chosen,'stable'); first=left(1);
 distinct=left(over(left)==over(first)&coverage(left)==coverage(first));
 distinct=distinct(arrayfun(@(j)all(cellfun(@(n)~isequal(n.visited,nodes{j}.visited),nodes(chosen))),distinct));
 if ~isempty(distinct), first=distinct(1); end
 if c.cost_guidance_enabled&&coverage(first)<35&&numel(chosen)==1
  % Retain a time-resource extreme beside the cheapest progressing prefix.
  % This is beam diversity, not a new objective or a dominance proof.
  peers=left(over(left)==over(first)&coverage(left)==coverage(first));
  advancing=peers(cellfun(@(n)n.last_gain>0,nodes(peers)));
  if ~isempty(advancing), peers=advancing; end
  times=cellfun(@(n)n.t,nodes(peers));
  [~,earliest]=sortrows([times(:),reshape(score(peers),[],1)]);
  early=peers(earliest(1));
  if nodes{early}.t<nodes{chosen(1)}.t-c.beam_diversity_scales(3)
   first=early;
  end
 end
 chosen(end+1)=first;
end
zero=find(cellfun(@(n)n.last_gain==0,nodes));
zero=setdiff(zero,chosen,'stable');
[~,zorder]=sort(over(zero)); zero=zero(zorder);
quota=max(0,min(count-numel(chosen),max(1,floor(width/4))));
chosen=[chosen zero(1:min(numel(zero),quota))];
while numel(chosen)<count
 left=setdiff(1:numel(nodes),chosen,'stable'); nearest=inf(size(left));
 if any(~over(left)), left=left(~over(left)); nearest=inf(size(left)); end
 for k=1:numel(left)
  nearest(k)=min(cellfun(@(n)distance(nodes{left(k)},n,c),nodes(chosen)));
 end
 [~,k]=max(nearest); chosen(end+1)=left(k); %#ok<AGROW>
end
beam=nodes(chosen);
end

function d=distance(a,b,c)
s=c.beam_diversity_scales;
d=max([norm(a.state(1:3)-b.state(1:3))/s(1), ...
 norm(a.state(4:6)-b.state(4:6))/s(2),abs(a.t-b.t)/s(3)]);
if ~isequal(a.visited,b.visited), d=max(d,2); end
end
