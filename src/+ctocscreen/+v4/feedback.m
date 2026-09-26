function [memory,weight,details]=feedback(memory,operation,node,c)
%FEEDBACK Bounded state-conditioned cost/plane observations, with deduplication.
if isempty(memory)
 memory=struct('keys',{{}},'negative',[],'positive',[],'observed',{{}},'positive_observed',{{}},'population',{{}}, ...
  'queries',0,'nonneutral',0,'updates',0,'positive_updates',0);
end
weight=1; details=struct('key','','loss',0);
if strcmp(operation,'evaporate'), memory.negative=(1-c.evaporation)*memory.negative; return; end
keys=features(node);
if strcmp(operation,'query')
 memory.queries=memory.queries+numel(keys); values=ones(numel(keys),1);
 for k=1:numel(keys)
  idx=find(strcmp(memory.keys,keys{k}),1);
  if ~isempty(idx)
   values(k)=max(.05,(1+memory.positive(idx))/(1+memory.negative(idx)));
   if abs(values(k)-1)>1e-12, memory.nonneutral=memory.nonneutral+1; end
  end
 end
 if ~isempty(values), weight=exp(mean(log(values))); end
 return
end
physical=ctocscreen.v4.controlKey(node.q);
seen=any(strcmp(memory.observed,physical));
eligible=node.actual.independent&&node.actual.passed&&node.actual.total_dv_km_s<=c.search_max_dv_km_s;
if seen&&(~eligible||any(strcmp(memory.positive_observed,physical))), return; end
if strcmp(node.actual.status,'propagation_failure')||~isfinite(node.actual.total_dv_km_s), return; end
if ~seen, memory.observed{end+1}=physical; end
J=node.actual.total_dv_km_s; costs=vecnorm(node.q.u,2,2);
loss=max(0,J/c.search_max_dv_km_s-1)*costs/max(J,eps);
if numel(node.actual.inclination_changes_deg)==numel(costs)
 loss=loss+(max(0,node.actual.inclination_changes_deg-c.plane_threshold_deg)/c.plane_threshold_deg).^2;
end
for k=1:numel(keys)
 idx=find(strcmp(memory.keys,keys{k}),1);
 if isempty(idx)
  memory.keys{end+1}=keys{k}; idx=numel(memory.keys); memory.negative(idx)=0; memory.positive(idx)=0;
 end
 if ~seen, memory.negative(idx)=min(c.negative_cap,memory.negative(idx)+min(c.negative_cap,loss(k))); end
end
memory.updates=memory.updates+1; details.loss=sum(loss);
if eligible
 memory.positive_observed{end+1}=physical;
 memory.population{end+1}=keys; memory.positive_updates=memory.positive_updates+1;
 if numel(memory.population)>c.positive_capacity, memory.population(1)=[]; end
 memory.positive(:)=0;
 for member=1:numel(memory.population)
  for k=1:numel(memory.population{member})
   idx=find(strcmp(memory.keys,memory.population{member}{k}),1);
   memory.positive(idx)=memory.positive(idx)+2/c.positive_capacity;
  end
 end
end
end
function keys=features(node)
keys=cell(numel(node.q.tau),1);
for k=1:numel(keys)
 x=ctocscreen.v4.stateAt(node.trace,node.q.tau(k),'pre');
 h=cross(x(1:3),x(4:6)); h=h/max(norm(h),eps);
 next=node.q.T; if k<numel(keys), next=node.q.tau(k+1); end
 seen=node.actual.witness_times_s<=node.q.tau(k)&node.actual.distance_km<=1;
 keys{k}=sprintf('%d,',[round(node.q.tau(k)/21600);round(norm(x(1:3))/5000); ...
  round(h*5);round(norm(x(4:6))/.5);round((next-node.q.tau(k))/7200); ...
  round(node.q.u(k,:).'/0.5);sum(seen(1:12));sum(seen(13:24));sum(seen(25:35))]);
end
end
