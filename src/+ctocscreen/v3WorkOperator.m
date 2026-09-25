function [state,operator]=v3WorkOperator(state,schedule,c)
%V3WORKOPERATOR Rotate work by lineage; a refined child cannot reset the cycle.
key=schedule.root_key;
index=find(strcmp(state.work_lineage_keys,key),1);
if isempty(index)
 state.work_lineage_keys{end+1}=key; index=numel(state.work_lineage_keys);
 state.work_lineage_attempts(index)=0;
end
phase=mod(state.work_lineage_attempts(index),3);
state.work_lineage_attempts(index)=state.work_lineage_attempts(index)+1;
operator='mutation';
if phase==0, operator='joint_refine';
elseif phase==1&&c.replan_enabled, operator='directed_replan'; end
end
