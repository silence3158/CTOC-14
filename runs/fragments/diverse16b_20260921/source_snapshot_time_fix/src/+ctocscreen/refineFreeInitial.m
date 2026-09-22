function out=refineFreeInitial(source,p,cfg)
%REFINEFREEINITIAL Joint six initial elements, all waits and encounter times.
% Fixed target order / Lambert families; no frozen initial-orbit components.
plan=source.plan; [s,r]=ctocscreen.rebuildSchedule(plan,p,cfg); assert(r.passed);
n=numel(plan.ids); assert(~any(plan.skip),'Reduced joint solve currently requires S1 topology.');
angleCenter=plan.initial_q(5:6);
offset=[p.re_km+600 zeros(1,5+2*n)]; scale=[10 .001 .001 pi 2*pi 2*pi repmat(10000,1,2*n)];
x0=([plan.initial_q plan.waits' plan.times']-offset)./scale;
lb=[-1 -.999999 -.999999 0 0 0 zeros(1,2*n)];
ub=[1 .999999 .999999 1 1 1 repmat(p.horizon_s/10000,1,2*n)];
lb(5:6)=angleCenter/(2*pi)-.5;ub(5:6)=angleCenter/(2*pi)+.5;
A=zeros(n,6+2*n); b=-60/10000*ones(n,1);
for k=1:n
 A(k,6+k)=1; A(k,6+n+k)=-1;
 if k>1, A(k,6+n+k-1)=1; end
end
out=struct('schedule',s,'evaluation',r,'plan',plan,'history',[0 r.total_dv_km_s], ...
 'evaluations',0,'invalid_evaluations',0,'failure_examples',{{}},'config',cfg,'source_hash',ctocscreen.implementationHash(), ...
 'scope','all six initial elements + all waits and arrival times; fixed order/branch families');
start=tic;
stageDeadline=.6*cfg.joint_seconds;
opt=optimoptions('fmincon','Algorithm','sqp','Display','off','MaxIterations',cfg.joint_iterations, ...
 'MaxFunctionEvaluations',cfg.joint_evaluations,'FiniteDifferenceStepSize',1e-6, ...
 'ConstraintTolerance',1e-9,'UseParallel',false,'OutputFcn',@stop);
% First couple the initial orbit to its first coast and encounter, then release
% all timing variables. This avoids a poor all-at-once direction hiding q0 gains.
active=[1:7 7+n];fixed=setdiff(1:numel(x0),active);
[~,~,out.block_exitflag,out.block_output]=fmincon(@blockObjective,x0(active), ...
 A(:,active),b-A(:,fixed)*x0(fixed)',[],[],lb(active),ub(active),@blockConstraints,opt);
encodedQ=out.plan.initial_q;
encodedQ(5:6)=angleCenter+mod(encodedQ(5:6)-angleCenter+pi,2*pi)-pi;
x0=([encodedQ out.plan.waits' out.plan.times']-offset)./scale;
stageDeadline=cfg.joint_seconds;
if toc(start)<cfg.joint_seconds
 [~,~,out.exitflag,out.output]=fmincon(@objective,x0,A,b,[],[],lb,ub,@constraints,opt);
else
 out.exitflag=out.block_exitflag;out.output=out.block_output;
end
out.elapsed_s=toc(start); out.initial_q_change=out.plan.initial_q-source.plan.initial_q;
out.independent=ctocscreen.propagateSchedule(out.schedule,p,true);
 function value=objective(x)
  out.evaluations=out.evaluations+1; v=x.*scale+offset;
  trial=plan; trial.initial_q=v(1:6); trial.waits=v(7:6+n)'; trial.times=v(7+n:end)';
  trial.initial_q(5:6)=mod(trial.initial_q(5:6),2*pi);
  [ss,rr]=ctocscreen.rebuildSchedule(trial,p,cfg); value=1e6;
  % Finite differences can cross the nonlinear eccentricity bound. The
  % numerical objective must remain smooth there; only valid rr enter archive.
  if isfield(rr,'build_completed')&&rr.build_completed
   value=sum(vecnorm(ss.delta_v_km_s,2,2));
  end
  if rr.passed
   value=rr.total_dv_km_s;
   if value<out.evaluation.total_dv_km_s
    out.plan=trial;out.schedule=ss;out.evaluation=rr;
    out.history(end+1,:)=[toc(start) value];
   end
  else
   out.invalid_evaluations=out.invalid_evaluations+1;
   if numel(out.failure_examples)<5,out.failure_examples{end+1}=rr;end
  end
 end
 function [c,eq]=constraints(x), c=x(2)^2+x(3)^2-.999999^2; eq=[];end
 function value=blockObjective(y),xx=x0;xx(active)=y;value=objective(xx);end
 function [c,eq]=blockConstraints(y),xx=x0;xx(active)=y;[c,eq]=constraints(xx);end
 function yes=stop(~,~,~), yes=toc(start)>=stageDeadline;end
end
