function beam=selectBeam(nodes,c)
%SELECTBEAM Keep actual states, coverage combinations, and root representatives.
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
count=cellfun(@(n)n.actual.visit_count,valid); cost=cellfun(@(n)n.actual.total_dv_km_s,valid);
[~,order]=sortrows([-count(:),cost(:)]); selected=order(1); root=cellfun(@(n)n.root_id,valid);
for k=order(:).'
 if numel(selected)>=min(c.beam_width,ceil(c.beam_width/2)), break; end
 if ~ismember(root(k),root(selected)), selected(end+1)=k; end %#ok<AGROW>
end
F=zeros(numel(valid),43);
for k=1:numel(valid)
 n=valid{k}; x=n.actual.final_state;
 F(k,:)=[n.q.T/43200,x(1:3).'/10000,x(4:6).'/0.5,2*(n.actual.distance_km<=1).',n.root_id/4];
end
while numel(selected)<min(c.beam_width,numel(valid))
 remaining=setdiff(1:numel(valid),selected,'stable'); score=zeros(numel(remaining),1);
 for j=1:numel(remaining)
  k=remaining(j); diversity=min(vecnorm(F(selected,:)-F(k,:),2,2));
  score(j)=count(k)-.05*cost(k)+.4*min(10,diversity);
 end
 [~,j]=max(score); selected(end+1)=remaining(j); %#ok<AGROW>
end
beam=valid(selected);
end
