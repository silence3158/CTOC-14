function f=refineFragment(xin,tin,ids,z0,limits,p,cfg)
%REFINEFRAGMENT S2 direct shooting: [wait,dv(3),dtA,dtAB], no impulse at A.
f=struct('passed',false,'status','not_run','evaluations',0,'exitflag',NaN);
scale=[10000 1 1 1 10000 10000]; start=tic;
lb=limits(1,:)./scale; ub=limits(2,:)./scale; y0=min(ub,max(lb,z0./scale));
opt=optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',cfg.fragment_evaluations, ...
 'MaxIterations',120,'FunctionTolerance',1e-12,'OptimalityTolerance',1e-10, ...
 'StepTolerance',1e-12,'FiniteDifferenceType','central','UseParallel',false,'OutputFcn',@stop);
try
 residual(y0); % Fail before solving on invalid initial data.
 [y,~,~,f.exitflag,f.output]=lsqnonlin(@residual,y0,lb,ub,opt);
 z=y.*scale; [rr,xa,xb,xplus]=residual(y);
 a=ctocscreen.checkArc(xin,z(1),p); b=ctocscreen.checkArc(xplus,z(5)+z(6),p);
 f.z=z; f.event_times_s=tin+z(1)+[z(5);z(5)+z(6)];
 f.event_distances_km=[norm(rr(1:3));norm(rr(4:6))]*1000;
 f.state_A=xa; f.state_out=xb; f.min_altitude_km=min(a.min_altitude_km,b.min_altitude_km);
 f.passed=all(f.event_distances_km<=cfg.fragment_tolerance_km)&& ...
  strcmp(a.status,'ok')&&strcmp(b.status,'ok')&&f.min_altitude_km>=200&& ...
  f.event_times_s(2)<=limits(3,1)&&f.event_times_s(2)<=p.horizon_s;
 f.status='infeasible'; if f.passed, f.status='fragment_screened'; end
catch err
 f.status='solver_failure'; f.failure_reason=err.message;
end
f.elapsed_s=toc(start);
 function [rr,xa,xb,xp]=residual(y)
  f.evaluations=f.evaluations+1; z=y.*scale;
  [xp,info]=ctocscreen.propagateTwoBody(xin,z(1),p.mu_km3_s2); assert(strcmp(info.status,'ok'));
  xp(4:6)=xp(4:6)+z(2:4);
  [xa,ia]=ctocscreen.propagateTwoBody(xp,z(5),p.mu_km3_s2);
  [xb,ib]=ctocscreen.propagateTwoBody(xa,z(6),p.mu_km3_s2);
  assert(strcmp(ia.status,'ok')&&strcmp(ib.status,'ok'));
  targets=ctocscreen.targetStates(p,ids,tin+z(1)+[z(5);z(5)+z(6)]);
  rr=[xa(1:3)-targets(1,1:3),xb(1:3)-targets(2,1:3)]'/1000;
 end
 function yes=stop(~,~,~), yes=toc(start)>=cfg.fragment_seconds; end
end
