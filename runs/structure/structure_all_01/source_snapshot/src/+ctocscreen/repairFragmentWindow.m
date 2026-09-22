function [bestPlan,bestS,bestR,report]=repairFragmentWindow(plan,s,r,p,cfg,stream)
%REPAIRFRAGMENTWINDOW S1 wait and encounter-time SQP, full suffix reevaluated.
bestPlan=plan; bestS=s; bestR=r; report=struct('evaluations',0,'exitflag',NaN);
n=numel(plan.ids); good=find(~plan.skip);
% Do not change a fixed S2 departure or its immediately preceding boundary.
good=good(all(isnan(plan.fixed_dv(good,:)),2));
good=good(good==n | ~ismember(good+1,find(any(isfinite(plan.fixed_dv),2))));
if isempty(good), return; end
k=good(randi(stream,numel(good))); report.k=k;
t0=0; if k>1, t0=plan.times(k-1); end
tend=p.horizon_s; if k<n, tend=plan.times(k+1)-60; end
span=tend-t0; scale=10000;
y0=[plan.waits(k),plan.times(k)-t0]/scale;
lo=[0 max(60,.55*(plan.times(k)-t0))]/scale;
hi=[min(.6*span,20000),min(span,1.45*(plan.times(k)-t0))]/scale;
start=tic;
opt=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
 'MaxFunctionEvaluations',cfg.local_evaluations,'MaxIterations',20, ...
 'FiniteDifferenceStepSize',1e-5,'UseParallel',false,'OutputFcn',@stop);
try
 [~,~,report.exitflag,report.output]=fmincon(@objective,y0,[1 -1],-60/scale,[],[],lo,hi,[],opt);
catch err
 report.failure_reason=err.message;
end
report.elapsed_s=toc(start);
 function value=objective(y)
  report.evaluations=report.evaluations+1;
  trial=plan; trial.waits(k)=y(1)*scale; trial.times(k)=t0+y(2)*scale;
  [ss,rr]=ctocscreen.rebuildSchedule(trial,p,cfg);
  value=1e6; % Solver-only infeasibility penalty, never an archive value.
  if rr.passed
   value=rr.total_dv_km_s;
   if value<bestR.total_dv_km_s, bestPlan=trial; bestS=ss; bestR=rr; end
  end
 end
 function yes=stop(~,~,~), yes=toc(start)>=cfg.local_seconds; end
end
