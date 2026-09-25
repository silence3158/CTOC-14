function [queue,added]=v3QueueCandidates(queue,items,capacity)
%V3QUEUECANDIDATES FIFO admission preserves unprocessed independent candidates.
added=0;
for k=1:numel(items)
 item=items{k}; key=itemKey(item);
 if any(cellfun(@(q)strcmp(itemKey(q),key),queue)), continue; end
 if numel(queue)>=capacity, break; end
 queue{end+1}=item; added=added+1;
end
end

function key=itemKey(item)
parent=[]; if isfield(item,'parent'), parent=item.parent; end
key=ctocscreen.v3RecoveryKey(item.schedule,item.target_ids,parent);
end
