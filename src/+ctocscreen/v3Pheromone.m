function [value,table,detail]=v3Pheromone(operation,table,key,amount,c)
%V3PHEROMONE Separate full-task rewards and bounded P-ACO prefix guidance.
value=c.pheromone_baseline;
detail=struct('full_increment',0,'prefix_increment',0,'negative_penalty',0);
switch operation
 case 'get'
  if isKey(table,key), value=table(key); end
  detail.full_increment=max(0,value-c.pheromone_baseline);
  guidance=['prefix|' key];
  if isKey(table,guidance), detail.prefix_increment=table(guidance); end
  value=min(c.pheromone_ceiling,value+detail.prefix_increment);
  negative=['negative|' key];
  if isKey(table,negative), detail.negative_penalty=table(negative); end
  value=max(c.pheromone_floor,value/(1+detail.negative_penalty));
 case 'evaporate'
  kk=keys(table);
  for j=1:numel(kk)
   if startsWith(kk{j},'prefix|'), continue; end
   if startsWith(kk{j},'observed|'), continue; end
   if startsWith(kk{j},'negative|')
    table(kk{j})=(1-c.evaporation)*table(kk{j}); continue;
   end
   table(kk{j})=c.pheromone_baseline+(1-c.evaporation)*max(0,table(kk{j})-c.pheromone_baseline);
  end
 case 'deposit'
  if isKey(table,key), value=table(key); end
  validateattributes(amount,{'double'},{'scalar','nonnegative','finite'});
  value=min(c.pheromone_ceiling,max(c.pheromone_baseline,value)+amount);
  table(key)=value;
 case 'penalize'
  validateattributes(amount,{'double'},{'scalar','nonnegative','finite'});
  negative=['negative|' key]; previous=0;
  if isKey(table,negative), previous=table(negative); end
  table(negative)=min(c.policy_penalty_ceiling,previous+amount);
  value=max(c.pheromone_floor,c.pheromone_baseline/(1+table(negative)));
 otherwise
  error('ctocscreen:v3:pheromone','Unknown operation.');
end
end
