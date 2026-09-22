function [s,r]=rebuildArcPlan(plan,p,cfg)
%REBUILDARCPLAN Lambert eliminates terminal position, never resets any state.
s=struct('schema_version','fragment_schedule_v2','dynamics_id','two_body', ...
 'initial_q',plan.initial_q,'event_target_ids',plan.ids(:),'event_times_s',plan.times(:), ...
 'duration_s',plan.times(end),'maneuver_times_s',[],'delta_v_km_s',[]);
r=struct('passed',false,'total_dv_km_s',Inf,'completed',false);
try
 assert(all(diff(plan.times)>0)&&plan.times(1)>0&&plan.times(end)<=p.horizon_s);
 x=ctocscreen.initialState(plan.initial_q,p.mu_km3_s2);time=0;last=0;minimum=Inf;
 states=zeros(numel(plan.ids),6);residual=zeros(numel(plan.ids),3);pre=zeros(numel(plan.counts),6);
 for m=1:numel(plan.counts)
  ix=last+(1:plan.counts(m));arrival=plan.times(ix(end));wait=plan.waits(m);
  assert(wait>=0&&time+wait<plan.times(ix(1)));
  a=ctocscreen.checkArc(x,wait,p);assert(strcmp(a.status,'ok'));minimum=min(minimum,a.min_altitude_km);
  x=ctocscreen.propagateTwoBody(x,wait,p.mu_km3_s2);time=time+wait;
  target=ctocscreen.targetStates(p,plan.ids(ix(end)),arrival);
  bs=ctocscreen.enumerateBranches(x(1:3),target(1:3),arrival-time,p.mu_km3_s2,cfg.branch_policy);
  assert(~isempty(bs));ref=plan.reference_v(m,:);
  if isfield(plan,'use_greedy_branches')&&plan.use_greedy_branches,ref=x(4:6);end
  [~,order]=sort(arrayfun(@(b)norm(b.v_depart-ref),bs));
  found=false;
  for b=order(:)'
   a=ctocscreen.checkArc([x(1:3) bs(b).v_depart],arrival-time,p);
   if strcmp(a.status,'ok')&&a.min_altitude_km>=200,found=true;break;end
  end
  assert(found);pre(m,:)=x;u=bs(b).v_depart-x(4:6);x(4:6)=x(4:6)+u;
  s.maneuver_times_s(m,1)=time;s.delta_v_km_s(m,:)=u;
  for j=ix
   dt=plan.times(j)-time;a=ctocscreen.checkArc(x,dt,p);assert(strcmp(a.status,'ok'));
   minimum=min(minimum,a.min_altitude_km);x=ctocscreen.propagateTwoBody(x,dt,p.mu_km3_s2);
   time=plan.times(j);states(j,:)=x;target=ctocscreen.targetStates(p,plan.ids(j),time);
   residual(j,:)=x(1:3)-target(1:3);
  end
  last=ix(end);
 end
 q=plan.initial_q;validQ=q(1)>=p.re_km+590&&q(1)<=p.re_km+610&&hypot(q(2),q(3))<.001&&q(4)>=0&&q(4)<=pi;
 distances=vecnorm(residual,2,2);
 r=struct('completed',true,'passed',validQ&&minimum>=200&&all(distances<=cfg.fragment_tolerance_km)&& ...
  numel(unique(plan.ids))==size(p.states0,1),'total_dv_km_s',sum(vecnorm(s.delta_v_km_s,2,2)), ...
  'position_residual_km',residual,'event_distances_km',distances,'event_states',states, ...
  'preburn_states',pre,'min_altitude_km',minimum);
catch err
 r.failure_reason=err.message;r.failure_stack=err.stack;
end
end
