function key=v3RecoveryKey(s,ids,parent)
%V3RECOVERYKEY Identity of a recovery task, distinct from physical feedback.
if nargin<3, parent=[]; end
ids=sort(ids(:));
plan=s.witness_times_s; explicit=isfield(s,'visit_plan_times_s');
if explicit, plan=s.visit_plan_times_s; end
key=[ctocscreen.v3ScheduleKey(s),'|targets=',mat2str(ids,17), ...
 '|explicit_plan=',mat2str(explicit),'|times=',mat2str(plan(ids),17)];
if ~isempty(parent)
 key=[key,'|parent=',ctocscreen.v3ScheduleKey(parent.schedule)];
end
end
