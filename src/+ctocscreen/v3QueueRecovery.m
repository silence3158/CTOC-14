function state=v3QueueRecovery(state,schedule,ids,parent,stalls,c)
%V3QUEUERECOVERY Separate unchanged recovery retries from directed rebuilding.
if isempty(ids), return; end
item=struct('schedule',schedule,'target_ids',ids,'parent',parent,'stalls',stalls);
if stalls>=c.recovery_stall_limit
 key=ctocscreen.v3RecoveryKey(schedule,ids,parent);
 if any(strcmp(state.replanned_keys,key)), return; end
 [state.replan_queue,~]=ctocscreen.v3QueueCandidates(state.replan_queue,{item},c.archive_size);
else
 [state.recovery_pool,~]=ctocscreen.v3QueueCandidates(state.recovery_pool,{item},c.archive_size);
end
end
