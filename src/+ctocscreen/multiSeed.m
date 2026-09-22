function [xin,tin,ids,z,finish]=multiSeed(s,r,k,counts,p,stream)
%MULTISEED Warm initialization only; returned trajectory need not be feasible.
N=sum(counts); ids=s.event_target_ids(k:k+N-1); tin=0;
xin=ctocscreen.initialState(s.initial_q,p.mu_km3_s2);
if k>1, tin=s.event_times_s(k-1); xin=r.event_states(k-1,:); end
finish=p.horizon_s; if k+N<=numel(s.event_times_s), finish=s.event_times_s(k+N)-60; end
visitTimes=s.event_times_s(k:k+N-1); eventTimes=[]; u=[]; visited=0; previous=tin;
for m=1:numel(counts)
 nextVisit=visitTimes(visited+1); burn=previous+.12*(nextVisit-previous);
 eventTimes(end+1)=burn; %#ok<AGROW>
 j=find(s.maneuver_times_s<nextVisit,1,'last');
 v=s.delta_v_km_s(j,:);
 if counts(m)==0, v=.15*v+.03*randn(stream,1,3); end
 u(end+1,:)=v; %#ok<AGROW>
 for a=1:counts(m), visited=visited+1; eventTimes(end+1)=visitTimes(visited); end %#ok<AGROW>
 previous=eventTimes(end);
end
z=[diff([tin eventTimes]) reshape(u',1,[])];
end
