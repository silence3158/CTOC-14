function [source,log]=kickArcSeed(base,p,cfg,stream)
%KICKARCSEED Nonlocal basin change; a feasible worse seed is deliberately allowed.
start=tic;source=[];log=struct('attempts',0,'screened',0,'independent_rejected',0,'mode',[]);
old=ctocscreen.arcPlan(base.schedule,base.evaluation);N=numel(old.ids);ends=cumsum(old.counts);
pl=old;pl.counts=ones(N,1);pl.waits=zeros(N,1);pl.reference_v=zeros(N,3);
for j=1:N
 m=find(ends>=j,1);pl.reference_v(j,:)=old.reference_v(m,:);
 if j==1||ends(m)-old.counts(m)+1==j,pl.waits(j)=old.waits(m);end
end
pool={};costs=[];
while log.attempts<cfg.kick_attempts&&toc(start)<cfg.kick_seconds
 log.attempts=log.attempts+1;trial=pl;mode=mod(log.attempts+randi(stream,3),3)+1;
 % Every proposal changes the sequence; temporal kicks are additional, not substitutes.
 first=randi(stream,[2 N-2]);last=min(N,first+randi(stream,[2 8]));
 switch mode
  case 1
   trial.ids([first last])=trial.ids([last first]);
  case 2
   block=first:min(N-1,first+randi(stream,[1 3]));remaining=setdiff(1:N,block,'stable');
   at=randi(stream,[2 numel(remaining)]);order=[remaining(1:at-1) block remaining(at:end)];trial.ids=pl.ids(order);
  case 3
   trial.ids(first:last)=flipud(trial.ids(first:last));
 end
 if isequal(trial.ids,pl.ids),continue;end
 if rand(stream)<cfg.kick_time_probability
  dt=diff([0;trial.times]);dt=dt.*exp(cfg.kick_time_sigma*randn(stream,N,1));
  total=min(p.horizon_s-60,pl.times(end)*(cfg.kick_duration_range(1)+diff(cfg.kick_duration_range)*rand(stream)));
  dt=60+(total-60*N)*dt/sum(dt);trial.times=cumsum(dt);trial.waits=min(trial.waits,dt*.6);
 end
 if rand(stream)<.3
  ix=randperm(stream,N,min(6,N));trial.reference_v(ix,:)=trial.reference_v(ix,:)+.6*randn(stream,numel(ix),3);
 end
 [s,r]=ctocscreen.rebuildArcPlan(trial,p,cfg);
 if r.passed&&r.total_dv_km_s<base.independent.total_dv_km_s*cfg.kick_cost_ratio
  log.screened=log.screened+1;pool{end+1}=struct('schedule',s,'evaluation',r,'arc_plan',trial,'mode',mode); %#ok<AGROW>
  costs(end+1)=r.total_dv_km_s; %#ok<AGROW>
 end
end
[~,order]=sort(costs);
for ii=order(:)'
 candidate=pool{ii};ir=ctocscreen.propagateSchedule(candidate.schedule,p,true);
 if ir.passed
  nr=ctocscreen.propagateSchedule(candidate.schedule,p,false);assert(nr.passed);
  candidate.evaluation=nr;candidate.independent=ir;source=candidate;log.mode=candidate.mode;break;
 end
 log.independent_rejected=log.independent_rejected+1;
end
log.elapsed_s=toc(start);
if ~isempty(source)
 log.order_hamming=sum(source.schedule.event_target_ids~=base.schedule.event_target_ids);
 log.time_rms_s=sqrt(mean((source.schedule.event_times_s-base.schedule.event_times_s).^2));
 log.starting_dv=source.independent.total_dv_km_s;
end
end
