function folders=run_v3_batch(label,budget_s,seeds,workers,options)
%RUN_V3_BATCH Independent seeds, same solver; explicit serial fallback.
if nargin<5, options=struct(); end
validateattributes(seeds,{'double'},{'vector','integer','nonnegative','finite'});
validateattributes(workers,{'double'},{'scalar','integer','positive'});
assert(numel(unique(seeds))==numel(seeds),'ctocscreen:v3:seeds','Seeds must be distinct.');
folders=cell(numel(seeds),1); parallel=false;
if workers>1
 try
  assert(license('test','Distrib_Computing_Toolbox'));
  pool=gcp('nocreate');
  if isempty(pool), parpool('local',workers); end
  parallel=true;
 catch err
  warning('ctocscreen:v3:serialFallback','Parallel unavailable; same algorithm runs serially: %s',err.message);
 end
end
if parallel
 parfor (k=1:numel(seeds),workers)
  folders{k}=run_v3_search(sprintf('%s_seed_%d',label,seeds(k)),budget_s,seeds(k),options);
 end
else
 for k=1:numel(seeds)
  folders{k}=run_v3_search(sprintf('%s_seed_%d',label,seeds(k)),budget_s,seeds(k),options);
 end
end
end
