function [table,population,report]=v3PrefixReinforce(table,population,node,eph,c)
%V3PREFIXREINFORCE Finite P-ACO memory of checked prefixes, not full-task reward.
report=struct('admitted',false,'reason','unchecked_prefix','visit_count',sum(node.visited), ...
 'full_task_reward',false,'key_count',0,'population_size',numel(population));
s=node.schedule;
if node.J>c.search_max_dv_km_s, report.reason='user_cost_rejected_no_positive_credit'; return; end
if ~isfield(node,'completion_check')||~node.completion_check.passed, return; end
checked=node.completion_check;
if ~isfield(checked,'schedule_key')||~strcmp(checked.schedule_key,ctocscreen.v3ScheduleKey(s)) ...
 ||~strcmp(checked.witness_key,mat2str(s.witness_times_s,17)), return; end
visits=sum(isfinite(s.witness_times_s));
if visits<2||visits>=35, report.reason='outside_prefix_scope'; return; end
physicalKey=ctocscreen.v3ScheduleKey(s); root=sprintf('%.17g/',s.initial_q);
if any(cellfun(@(p)strcmp(p.physical_key,physicalKey),population))
 report.reason='duplicate_physical_prefix'; return;
end
family=find(cellfun(@(p)strcmp(p.root,root),population),1);
if ~isempty(family)
 old=population{family};
 if visits<old.visit_count||(visits==old.visit_count&&node.J>=old.J)
  report.reason='no_coverage_or_cost_progress'; return;
 end
end
if isfield(checked,'experience_keys')&&isfield(checked,'experience_config')&& ...
 isequal(checked.experience_config.pheromone_state_scales,c.pheromone_state_scales)
 experience=checked.experience_keys;
else
 experience=ctocscreen.v3ScheduleExperience(s,eph,c,'prefix_witnesses');
end
if isempty(experience), report.reason='no_state_actions'; return; end
item=struct('physical_key',physicalKey,'root',root,'initial_q',s.initial_q, ...
 'visit_count',visits,'J',node.J,'keys',{experience});
if ~isempty(family), population(family)=[]; end
population{end+1}=item;
if numel(population)>c.prefix_memory_size, population(1)=[]; end
% Membership defines the guidance exactly. Eviction removes every old contribution.
oldKeys=keys(table); oldKeys=oldKeys(startsWith(oldKeys,'prefix|'));
if ~isempty(oldKeys), remove(table,oldKeys); end
increment=c.prefix_pheromone_gain/c.prefix_memory_size;
for k=1:numel(population)
 for j=1:numel(population{k}.keys)
  key=['prefix|' population{k}.keys{j}]; value=0;
  if isKey(table,key), value=table(key); end
  table(key)=value+increment;
 end
end
report.admitted=true; report.reason='population_prefix_guidance';
report.visit_count=visits; report.key_count=numel(experience); report.population_size=numel(population);
end
