function [cands,report]=crossing(node,eph,c,ev,info,stream,budget)
%CROSSING Timed Lambert proposals to target plane-crossing events (A2).
% Departure times: node end, apses and uniform waits along the ballistic
% continuation (a stored solution, so no burn is implied by a wait). Each
% proposal is a two-body Lambert seed; J2 correction and replay decide reality.
clock=tic; m=eph.model; cands=struct('target',{},'departure',{},'arrival',{},'x',{},'v',{},'dv',{},'score',{},'kind',{});
report=struct('events',numel(ev),'lambert',0,'branches',0,'seconds',0);
if isempty(ev)||isempty(info.sol), return; end
T=node.q.T; sol=info.sol;
% Time is a shared resource: arrivals later than the remaining pace
% (remaining time / remaining targets) pay a price in km/s (search heuristic).
% The pace keeps a reserve so the last targets do not inherit an exhausted
% horizon; late in the mission the price rises as the reserve is consumed.
left=sum(node.actual.distance_km>1);
pace=max(c.min_flight_s,(1-c.time_reserve)*(m.horizon_s-T)/max(1,left));
timeCost=@(arrival)c.time_price_km_s*max(0,(arrival-T)/pace-1);
% Rank events by the cheap estimate; keep the best few per target for spread.
[~,order]=sort([ev.dv_est]+timeCost([ev.time])); picked=[]; perTarget=zeros(35,1);
for k=order
 if perTarget(ev(k).id)>=c.events_per_target, continue; end
 picked(end+1)=k; perTarget(ev(k).id)=perTarget(ev(k).id)+1; %#ok<AGROW>
 if numel(picked)>=c.event_candidates, break; end
end
% Candidate departures: now, apses of the continuation, and uniform waits.
tg=linspace(T,info.t_end,max(3,ceil((info.t_end-T)/300)));
X=deval(sol,tg); rr=vecnorm(X(1:3,:)); rdot=sum(X(1:3,:).*X(4:6,:),1)./rr;
k=find(rdot(1:end-1).*rdot(2:end)<=0); apses=tg(k)+(tg(k+1)-tg(k)).*rdot(k)./(rdot(k)-rdot(k+1));
policy=struct('max_revolutions',c.max_revolutions,'endpoint_tol_km',.005);
for k=picked
 e=ev(k); latest=e.time-c.min_flight_s;
 departures=[T,apses(apses<latest),T+(latest-T)*[.25 .5]];
 departures=unique(departures(departures>=T&departures<=latest));
 if numel(departures)>c.departures_per_event
  % Keep now, the latest apses and a waited point.
  departures=departures(unique(round(linspace(1,numel(departures),c.departures_per_event))));
 end
 rj=ctocscreen.v3QueryTargets(eph,e.id,e.time);
 for tau=departures
  if toc(clock)>=budget, break; end
  if tau==T, xp=node.actual.final_state; else, xp=deval(sol,tau); xp=xp(1:6); end
  report.lambert=report.lambert+1;
  try
   branches=ctocscreen.v3LambertBranches(xp(1:3),rj,e.time-tau,m.mu,policy);
  catch
   continue
  end
  for b=1:numel(branches)
   v=branches(b).v_depart(:); dv=norm(v-xp(4:6));
   h1=cross(xp(1:3),xp(4:6)); h2=cross(xp(1:3),v);
   inc=abs(acosd(max(-1,min(1,h1(3)/norm(h1))))-acosd(max(-1,min(1,h2(3)/norm(h2)))));
   penalty=(max(0,inc-c.plane_threshold_deg)/c.plane_threshold_deg)^2;
   cands(end+1)=struct('target',e.id,'departure',tau,'arrival',e.time,'x',xp,'v',v, ...
    'dv',dv,'score',dv+c.plane_weight*min(4,penalty)+timeCost(e.time),'kind','crossing'); %#ok<AGROW>
   report.branches=report.branches+1;
  end
 end
 if toc(clock)>=budget, break; end
end
report.seconds=toc(clock);
end
