function [f,log]=refineMultiFragment(xin,tin,ids,counts,z0,p,cfg,finish)
%REFINEMULTIFRAGMENT Generic M-burn/N-target direct shooting and constrained polish.
L=numel(counts)+numel(ids); scale=[repmat(10000,1,L) ones(1,3*numel(counts))];
span=finish-tin; lb=[0 repmat(30,1,L-1) z0(L+1:end)-cfg.multi_velocity_window];
ub=[repmat(span,1,L) z0(L+1:end)+cfg.multi_velocity_window];
% Bound each interval around its seed, without mistaking it for a mission bound.
ub(1:L)=min(ub(1:L),max(600,z0(1:L)*2+3000));
start=tic; log=struct('counts',counts,'evaluations',0,'exitflag',NaN,'initial_cost',Inf);
f=struct('passed',false,'status','solver_failure');
try
 algorithm='trust-region-reflective';
 if 3*numel(ids)<numel(z0), algorithm='levenberg-marquardt';end
 opt=optimoptions('lsqnonlin','Algorithm',algorithm,'Display','off','MaxFunctionEvaluations',cfg.multi_evaluations, ...
  'MaxIterations',180,'FunctionTolerance',1e-12,'OptimalityTolerance',1e-10, ...
  'StepTolerance',1e-11,'FiniteDifferenceType','central','UseParallel',false,'OutputFcn',@stop);
 y0=min(ub,max(lb,z0))./scale; residual(y0);
 [y,~,~,log.exitflag,log.output]=lsqnonlin(@residual,y0,lb./scale,ub./scale,opt);
 f=ctocscreen.replayMultiFragment(xin,tin,ids,counts,y.*scale,p);
 f.passed=f.passed&&f.t_out_s<=finish;
 if ~f.passed,f.status='fragment_infeasible';end
 log.initial_cost=f.total_dv_km_s;
 if f.passed&&toc(start)<cfg.multi_seconds*.8
  % Feasible starting point: optimize actual burn cost under explicit constraints.
  best=f; bestY=y;
  op=optimoptions('fmincon','Algorithm','sqp','Display','off','MaxIterations',20, ...
   'MaxFunctionEvaluations',cfg.multi_polish_evaluations,'UseParallel',false, ...
   'ConstraintTolerance',1e-9,'OutputFcn',@stop);
  [~,~,log.polish_exitflag,log.polish_output]=fmincon(@cost,y, ...
   [scale(1:L) zeros(1,numel(y)-L)],span,[],[],lb./scale,ub./scale,@constraints,op);
  f=best; y=bestY;
 end
 f.z=y.*scale;
catch err
 log.failure_reason=err.message;
end
log.elapsed_s=toc(start);
 function rr=residual(y)
  log.evaluations=log.evaluations+1;
  ff=ctocscreen.replayMultiFragment(xin,tin,ids,counts,y.*scale,p);
  % Position-only feasibility search. Height/time are explicitly checked afterwards.
  rr=reshape(ff.position_residual_km',[],1)/1000;
 end
 function v=cost(y)
  log.evaluations=log.evaluations+1;
  ff=ctocscreen.replayMultiFragment(xin,tin,ids,counts,y.*scale,p); v=ff.total_dv_km_s;
  if ff.passed&&ff.t_out_s<=finish&&v<best.total_dv_km_s, best=ff;bestY=y;end
 end
 function [c,ceq]=constraints(y)
  ff=ctocscreen.replayMultiFragment(xin,tin,ids,counts,y.*scale,p);
  c=[(ff.event_distances_km-.02)/1000;(200-ff.min_altitude_km)/1000]; ceq=[];
 end
 function yes=stop(~,~,~), yes=toc(start)>=cfg.multi_seconds;end
end
