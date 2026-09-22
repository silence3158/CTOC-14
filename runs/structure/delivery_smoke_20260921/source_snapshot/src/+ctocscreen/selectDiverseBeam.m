function selected=selectDiverseBeam(nodes,width)
%SELECTDIVERSEBEAM Heuristic diversity buckets; no claimed dominance proof.
selected={}; if isempty(nodes), return; end
cost=cellfun(@(a)a.cost,nodes); [~,ix]=sort(cost); keys={}; deferred={};
for j=ix
 a=nodes{j}; h=cross(a.state(1:3),a.state(4:6)); h=h/max(norm(h),eps);
 key=sprintf('%s_%d_%d_%d_%d',sprintf('%d',a.mask),round(a.time/10000),round(h(1)*4),round(h(2)*4),round(h(3)*4));
 if ~ismember(key,keys)
  selected{end+1}=a; keys{end+1}=key; %#ok<AGROW>
 else
  deferred{end+1}=a; %#ok<AGROW>
 end
 if numel(selected)>=width, return;end
end
for j=1:min(width-numel(selected),numel(deferred)), selected{end+1}=deferred{j};end
end
