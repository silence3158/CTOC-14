function source=v3ChooseSource(attempted,hasWork,hasRecovery)
%V3CHOOSESOURCE Rotate attempts, including failed and skipped source slots.
source='fresh';
if hasWork&&mod(attempted,2)==1, source='work';
elseif hasRecovery&&mod(attempted,3)==2, source='recovery'; end
if ~hasWork&&hasRecovery&&mod(attempted,2)==1, source='recovery'; end
end
