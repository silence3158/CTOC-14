function [table,ledger,report]=v3Reinforce(table,ledger,features,keys,J,c)
%V3REINFORCE Verified full-task reward; near copies earn only incremental gain.
validateattributes(J,{'double'},{'scalar','nonnegative','finite'});
dist=cellfun(@(entry)ctocscreen.v3FeatureDistance(entry.features,features,c),ledger);
close=find(dist<1); quality=c.pheromone_cost_scale_km_s/(c.pheromone_cost_scale_km_s+J);
if isfield(features,'family')
 same=find(cellfun(@(entry)isfield(entry.features,'family')&&strcmp(entry.features.family,features.family),ledger));
 close=unique([close same]);
end
amount=c.pheromone_gain*quality; reason='new_trajectory_family';
if ~isempty(close)
 [~,best]=min(cellfun(@(entry)entry.best_J,ledger(close))); index=close(best);
 old=ledger{index};
 if J>=old.best_J
  amount=0; reason='no_new_family_or_cost_gain';
 else
  oldQuality=c.pheromone_cost_scale_km_s/(c.pheromone_cost_scale_km_s+old.best_J);
  amount=c.pheromone_gain*(quality-oldQuality); reason='incremental_cost_gain';
  ledger{index}.best_J=J;
 end
else
 ledger{end+1}=struct('features',features,'best_J',J);
end
keys=unique(keys);
if amount>0
 for k=1:numel(keys)
  [~,table]=ctocscreen.v3Pheromone('deposit',table,keys{k},amount,c);
 end
end
report=struct('amount_per_key',amount,'key_count',numel(keys),'reason',reason,'full_task_dv_km_s',J);
end
