function [plans,log]=lowSpeedSeeds(source,base,core,p,cfg)
%LOWSPEEDSEEDS Guided P-topology starts. Discovery time belongs to L's budget.
start=tic;plans={};rows=zeros(0,7);scores=[];
m=core.first_arc;next=m+1;
burn=source.schedule.maneuver_times_s(m);finish=base.times(core.first_event);
x=source.evaluation.preburn_states(m,:);u=base.zero_delta_v(m,:);
R=x(1:3)/norm(x(1:3));H=cross(x(1:3),x(4:6));H=H/norm(H);T=cross(H,R);
directions=[zeros(1,3);.25*R;-.25*R;.25*T;-.25*T;.25*H;-.25*H;.6*T;-.6*T;.6*H;-.6*H];
baseline=ctocscreen.propagateTwoBody([x(1:3) x(4:6)+u],base.waits(next),p.mu_km3_s2);
vref=norm(baseline(4:6));costref=source.independent.total_dv_km_s;
% Later and interior points first; all lie before the first core encounter.
fractions=[.7 .4 .9 .2 .55 .8];attempts=0;failures=0;
for f=fractions
 for k=1:size(directions,1)
  if toc(start)>=cfg.low_speed_seconds||isfile(cfg.stop_file),break;end
  attempts=attempts+1;pl=base;pl.waits(next)=max(cfg.minimum_gap_s,f*(finish-burn));
  pl.zero_delta_v(m,:)=u+directions(k,:);
  if isfield(pl,'locked_branch_ids'),pl=rmfield(pl,'locked_branch_ids');end
  [s,r]=ctocscreen.rebuildArcPlan(pl,p,cfg);
  if ~r.passed,failures=failures+1;continue;end
  before=r.preburn_states(next,4:6);after=before+s.delta_v_km_s(next,:);
  speed=norm(before);angle=acos(max(-1,min(1,dot(before,after)/(norm(before)*norm(after)))));
  % Actual reconnection and whole-tour cost screen geometry; speed alone is
  % not a certificate of a useful turn. The ratio cap is a search heuristic.
  if speed>=.98*vref||r.total_dv_km_s>cfg.low_speed_cost_ratio*costref,continue;end
  score=speed/vref+r.total_dv_km_s/costref+norm(s.delta_v_km_s(next,:))/costref;
  pl.locked_branch_ids=r.branch_ids;plans{end+1}=pl;scores(end+1)=score; %#ok<AGROW>
  rows(end+1,:)=[pl.waits(next) speed angle*180/pi norm(s.delta_v_km_s(next,:)) r.total_dv_km_s score k]; %#ok<AGROW>
 end
 if toc(start)>=cfg.low_speed_seconds||isfile(cfg.stop_file),break;end
end
[~,order]=sort(scores);order=order(1:min(cfg.low_speed_keep,numel(order)));
plans=plans(order);rows=rows(order,:);
log=struct('elapsed_s',toc(start),'attempts',attempts,'infeasible',failures, ...
 'reference_speed_km_s',vref,'retained',numel(plans),'candidates',rows, ...
 'columns',{{'coast_s','speed_km_s','turn_deg','reconnection_dv_km_s','total_dv_km_s','score','direction_index'}}, ...
 'scope','nominal seeds only; final acceptance requires independent whole-mission replay');
end
