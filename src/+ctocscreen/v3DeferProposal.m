function state=v3DeferProposal(state,source,item)
%V3DEFERPROPOSAL Preserve the active proposal when its joint call cannot start.
switch source
 case 'initial'
  state.pending_initial=[{item.schedule},state.pending_initial];
 case 'fresh'
  state.fresh_queue=[{item},state.fresh_queue];
 otherwise
  % One reserved active slot may exceed the normal admission capacity by one.
  % Deferral is not a failed recovery and must retain its parent and stalls.
  state.recovery_pool=[{item},state.recovery_pool];
end
end
