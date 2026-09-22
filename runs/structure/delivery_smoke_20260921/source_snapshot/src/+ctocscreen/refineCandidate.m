function [best,bestEv,report]=refineCandidate(c,p,cfg,remainingSeconds)
%REFINECANDIDATE Fixed-order/branch SQP over all 42 continuous variables.
if nargin<4, remainingSeconds=Inf; end
best=c; bestEv=ctocscreen.evaluate(c,p,cfg,true);
report=struct('exitflag',NaN,'output',struct(),'evaluations',0);
assert(license('test','Optimization_Toolbox'),'Optimization Toolbox required for configured SQP.');
scale=[10 .001 .001 pi 2*pi 2*pi repmat(86400,1,36)];
offset=[p.re_km+600 zeros(1,41)];
x0=([c.initial_q c.wait_s c.tof_s]-offset)./scale;
lb=([-10 -.000999 -.000999 0 0 0 0 repmat(cfg.min_tof_s,1,35)])./scale;
ub=([10 .000999 .000999 pi 2*pi 2*pi repmat(p.horizon_s,1,36)])./scale;
A=[zeros(1,6) scale(7:42)]; start=tic;
opt=optimoptions('fmincon','Algorithm','sqp','Display','none', ...
 'MaxFunctionEvaluations',cfg.refine_evaluations,'MaxIterations',cfg.refine_iterations, ...
 'FiniteDifferenceStepSize',1e-5,'ConstraintTolerance',1e-9,'UseParallel',false, ...
 'OutputFcn',@stop);
[~,~,report.exitflag,report.output]=fmincon(@objective,x0,A,p.horizon_s,[],[],lb,ub,@constraints,opt);
 function value=objective(x)
  trial=decode(x); ev=ctocscreen.evaluate(trial,p,cfg,true);
  report.evaluations=report.evaluations+1;
  if strcmp(ev.status,'two_body_verified')
   value=ev.total_dv_km_s;
   if value<bestEv.total_dv_km_s, best=trial; bestEv=ev; end
  else
   value=1e6; % Internal solver penalty only, never archived as feasible.
  end
 end
 function trial=decode(x)
  v=x.*scale+offset; trial=c; trial.initial_q=v(1:6); trial.wait_s=v(7); trial.tof_s=v(8:42);
 end
 function [ineq,eq]=constraints(x)
  ineq=x(2)^2+x(3)^2-.999999^2; eq=[];
 end
 function yes=stop(~,~,~)
  yes=toc(start)>=remainingSeconds;
 end
end
