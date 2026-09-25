function [pool,report]=v3SelectWorkPool(pool,item,c,workingMove)
%V3SELECTWORKPOOL Called only after independent verification and SA decision.
% The accepted SA state survives; historical best belongs to the export archive.
if nargin<4, workingMove=false; end
if isfield(item.features,'family')
 family=find(cellfun(@(w)isfield(w.features,'family')&&strcmp(w.features.family,item.features.family),pool));
 if numel(family)>=c.work_family_quota
  [worst,index]=max(cellfun(@(w)w.J,pool(family)));
  if ~workingMove&&item.J>=worst
   report=struct('retained',false,'near_replaced',0,'minimum_distance',Inf,'reason','lineage_structure_quota'); return;
  end
  pool(family(index))=[];
 end
end
distance=cellfun(@(w)ctocscreen.v3FeatureDistance(w.features,item.features,c),pool);
close=find(distance<1);
report=struct('retained',false,'near_replaced',0,'minimum_distance',Inf,'reason','');
if ~isempty(distance), report.minimum_distance=min(distance); end
if ~isempty(close)&&~workingMove&&min(cellfun(@(w)w.J,pool(close)))<=item.J
 report.reason='near_copy_without_improvement'; return;
end
report.near_replaced=numel(close);
pool(close)=[]; pool{end+1}=item; incoming=numel(pool);
if numel(pool)>c.archive_size
 costs=cellfun(@(w)w.J,pool); [~,costOrder]=sort(costs);
 if workingMove, chosen=incoming; else, chosen=costOrder(1); end
 if numel(chosen)<c.archive_size, chosen=unique([chosen costOrder(1)],'stable'); end
 while numel(chosen)<c.archive_size
  left=setdiff(costOrder,chosen,'stable'); nearest=inf(size(left));
  for k=1:numel(left)
   nearest(k)=min(cellfun(@(w)ctocscreen.v3FeatureDistance(pool{left(k)}.features,w.features,c),pool(chosen)));
  end
  [~,k]=max(nearest); chosen(end+1)=left(k); %#ok<AGROW>
 end
 report.retained=ismember(incoming,chosen); pool=pool(chosen);
else
 report.retained=true;
end
if report.retained, report.reason='retained_diverse_work'; else, report.reason='capacity_diversity_selection'; end
end
