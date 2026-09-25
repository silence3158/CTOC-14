function [keys,features,actions]=v3ScheduleExperience(s,eph,c,scope)
%V3SCHEDULEEXPERIENCE Physical state keys; verification is the caller's duty.
if nargin<4, scope='full'; end
[r,trace]=ctocscreen.v3Replay(s,eph,c,false,scope); keys={}; actions={};
assert(~strcmp(r.status,'propagation_failure'),'ctocscreen:v3:experience','Cannot credit failed propagation.');
if nargout>1, features=ctocscreen.v3ScheduleFeatures(s,eph,c,trace); end
[keys,actions]=ctocscreen.v3TraceExperience(s,trace,c);
end
