function report=verifyIndependent(candidate,problem,config)
%VERIFYINDEPENDENT Replay stored impulses with independent ode113 integration.
ev=ctocscreen.evaluate(candidate,problem,config,true);
report=struct('passed',false,'method','ode113_two_body','status','not_run', ...
 'failure_reason','','endpoint_errors_km',nan(1,35),'max_endpoint_error_km',Inf, ...
 'min_altitude_km',NaN,'total_dv_km_s',Inf,'unique_visit_count',0, ...
 'relative_tolerance',3e-14,'absolute_tolerance',1e-14);
if ~strcmp(ev.status,'two_body_verified')
 report.status=ev.status; report.failure_reason=ev.failure_reason; return;
end
try
mu=problem.mu_km3_s2; opt=odeset('RelTol',3e-14,'AbsTol',1e-14);
state=integrate(ev.initial_state,candidate.wait_s); errors=zeros(1,35);
minAltitude=Inf;
a=ctocscreen.checkArc(ev.initial_state,candidate.wait_s,problem);
minAltitude=min(minAltitude,a.min_altitude_km);
for k=1:35
 state(4:6)=state(4:6)+ev.delta_v(k,:);
 a=ctocscreen.checkArc(state,candidate.tof_s(k),problem);
 assert(strcmp(a.status,'ok'),'Arc check failed during independent replay.');
 minAltitude=min(minAltitude,a.min_altitude_km);
 state=integrate(state,candidate.tof_s(k));
 target=integrate(problem.states0(candidate.order(k),:),ev.arrive_times_s(k));
 errors(k)=norm(state(1:3)-target(1:3));
end
report=struct('passed',all(errors<=1)&&minAltitude>=200,'method','ode113_two_body', ...
 'status','completed','failure_reason','', ...
 'endpoint_errors_km',errors,'max_endpoint_error_km',max(errors), ...
 'min_altitude_km',minAltitude,'total_dv_km_s',sum(vecnorm(ev.delta_v,2,2)), ...
 'unique_visit_count',numel(unique(candidate.order(errors<=1))), ...
 'relative_tolerance',3e-14,'absolute_tolerance',1e-14);
catch err
 report.status='solver_failure'; report.failure_reason=err.message;
end
 function out=integrate(in,dt)
  if dt==0, out=in; return; end
  [times,sol]=ode113(@(~,s)[s(4:6);-mu*s(1:3)/norm(s(1:3))^3],[0 dt],in,opt);
  assert(times(end)==dt && all(isfinite(sol(end,:))), ...
   'ctocscreen:independent:incomplete','Independent integration did not reach the endpoint.');
  out=sol(end,:);
 end
end
