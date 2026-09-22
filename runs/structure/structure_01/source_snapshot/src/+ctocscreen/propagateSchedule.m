function r=propagateSchedule(s,p,independent)
%PROPAGATESCHEDULE Fixed impulses, actual states, all-target visit witnesses.
% Witness checks certify positive visits, not absence of additional encounters.
if nargin<3, independent=false; end
r=struct('passed',false,'status','invalid_input','failure_reason','', ...
 'total_dv_km_s',Inf,'unique_visit_count',0,'min_altitude_km',Inf, ...
 'event_distances_km',[],'event_states',[],'validation_level','proposal');
try
 q=s.initial_q; times=s.maneuver_times_s(:); dv=s.delta_v_km_s;
 et=s.event_times_s(:); ids=s.event_target_ids(:); T=s.duration_s;
 assert(strcmp(s.schema_version,'fragment_schedule_v2')&&strcmp(s.dynamics_id,'two_body'));
 assert(numel(q)==6&&all(isfinite(q))&&q(1)>=p.re_km+590&&q(1)<=p.re_km+610);
 assert(hypot(q(2),q(3))<.001&&q(4)>=0&&q(4)<=pi);
 assert(isfinite(T)&&T>0&&T<=p.horizon_s);
 assert(size(dv,1)==numel(times)&&size(dv,2)==3&&all(isfinite(dv(:))));
 assert(all(isfinite(times))&&all(times>=0)&all(times<T)&&all(diff(times)>0));
 assert(numel(et)==numel(ids)&&all(isfinite(et))&&all(et>=0)&all(et<=T));
 assert(all(ids==floor(ids))&&all(ids>=1)&all(ids<=size(p.states0,1)));
 knots=unique([0;times;et;T]); state=ctocscreen.initialState(q,p.mu_km3_s2);
 r.event_states=nan(numel(et),6); r.event_distances_km=inf(numel(et),1);
 r.preburn_states=nan(numel(times),6); r.arc_states=zeros(numel(knots)-1,6);
 r.arc_times_s=[knots(1:end-1) knots(2:end)];
 targets=p.states0; last=0;
 opt=odeset('RelTol',3e-14,'AbsTol',1e-14);
 for k=1:numel(knots)
  t=knots(k); dt=t-last;
  a=ctocscreen.checkArc(state,dt,p); assert(strcmp(a.status,'ok'),'Arc propagation failed.');
  r.min_altitude_km=min(r.min_altitude_km,a.min_altitude_km);
  state=advance(state,dt);
  if independent&&dt>0
   % Independent target integration shares the event time mesh, no resets.
   [~,yy]=ode113(@targetRhs,[0 dt],reshape(targets',[],1),opt);
   targets=reshape(yy(end,:),6,[])';
  end
  ii=find(et==t);
  if ~isempty(ii)
   if independent, ts=targets(ids(ii),:); else, ts=ctocscreen.targetStates(p,ids(ii),t); end
   r.event_states(ii,:)=repmat(state,numel(ii),1);
   r.event_distances_km(ii)=vecnorm(ts(:,1:3)-state(1:3),2,2);
  end
  b=find(times==t);
  if ~isempty(b), r.preburn_states(b,:)=state; state(4:6)=state(4:6)+dv(b,:); end
  if k<numel(knots), r.arc_states(k,:)=state; end
  last=t;
 end
 assert(all(isfinite(r.event_distances_km))&&all(isfinite(state)));
 r.unique_visit_count=numel(unique(ids(r.event_distances_km<=1)));
 r.total_dv_km_s=sum(vecnorm(dv,2,2)); r.final_state=state;
 r.passed=r.unique_visit_count==size(p.states0,1)&&r.min_altitude_km>=200;
 r.status='constraint_violation';
 if r.passed
  r.status='completed'; r.validation_level='trajectory_screened';
  if independent, r.validation_level='two_body_independent_verified'; end
 end
 r.independent=independent; r.relative_tolerance=3e-14; r.absolute_tolerance=1e-14;
 r.event_check='all-target positive witnesses; no claim of exhaustive first-entry times';
catch err
 r.status='solver_failure'; r.failure_reason=err.message; r.failure_stack=err.stack; r.passed=false; r.total_dv_km_s=Inf;
end
 function out=advance(in,dt)
  if dt==0, out=in; return; end
  if independent
   [tt,yy]=ode113(@(~,x)[x(4:6);-p.mu_km3_s2*x(1:3)/norm(x(1:3))^3],[0 dt],in,opt);
   assert(tt(end)==dt&&all(isfinite(yy(end,:)))); out=yy(end,:);
  else
   [out,info]=ctocscreen.propagateTwoBody(in,dt,p.mu_km3_s2); assert(strcmp(info.status,'ok'));
  end
 end
 function dy=targetRhs(~,y)
  x=reshape(y,6,[]); a=-p.mu_km3_s2*x(1:3,:)./vecnorm(x(1:3,:),2,1).^3;
  dy=reshape([x(4:6,:);a],[],1);
 end
end
