function out=refineArcWindow(plan,focus,p,cfg)
%REFINEARCWINDOW Joint incoming/core/outgoing times; full remaining-mission cost.
% Solver penalties are search-only. Only explicit physical replay may pass.
start=tic;N=numel(plan.ids);M=numel(plan.counts);ends=cumsum(plan.counts);begins=[1;ends(1:end-1)+1];
arcs=max(1,focus-1):min(M,focus+cfg.lookahead_arcs);
multi=find(plan.counts>1);arcs=unique([arcs(:);multi(:);max(1,multi(:)-1);min(M,multi(:)+1)]);
events=[];for a=arcs',events=[events begins(a):ends(a)];end %#ok<AGROW>
events=unique(events);active=[events(:);N+arcs(:)];base=[plan.times;plan.waits];y0=base(active)/10000;
gap=diff([0;plan.times;p.horizon_s+3600]);radius=cfg.time_radius*min(gap(1:N),gap(2:N+1));
lower=[max(1,plan.times-radius);zeros(M,1)];upper=[min(p.horizon_s,plan.times+radius);zeros(M,1)];
for m=1:M
 prev=0;if m>1,prev=plan.times(ends(m-1));end
 upper(N+m)=max(plan.waits(m)+600,plan.times(begins(m))-prev-30);
end
intermediate=setdiff(1:N,ends);cacheY=[];cacheR=[];cacheS=[];evals=0;
out=struct('plan',plan,'schedule',[],'evaluation',struct('passed',false),'evaluations',0, ...
 'initial_cost',Inf,'exitflag',NaN,'feasibility_exitflag',NaN,'active_events',events,'active_arcs',arcs);
bestCost=Inf;
try
 evaluate(y0);if cacheR.completed,out.initial_cost=cacheR.total_dv_km_s;end
 if ~isfield(cfg,'visit_radius_km')&&~isempty(intermediate)&&(~cacheR.completed||max(cacheR.event_distances_km)>.03)
  op=optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt', ...
   'MaxFunctionEvaluations',cfg.window_evaluations,'MaxIterations',100,'FunctionTolerance',1e-12, ...
   'StepTolerance',1e-10,'FiniteDifferenceType','central','UseParallel',false,'OutputFcn',@stopFeas);
  [y0,~,~,out.feasibility_exitflag]=lsqnonlin(@residual,y0,lower(active)/10000,upper(active)/10000,op);
 end
 if toc(start)<cfg.window_seconds
  op=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxFunctionEvaluations',cfg.window_evaluations, ...
   'MaxIterations',80,'FiniteDifferenceStepSize',1e-6,'ConstraintTolerance',1e-9, ...
   'StepTolerance',1e-9,'OptimalityTolerance',1e-6,'UseParallel',false,'OutputFcn',@stop);
  [~,~,out.exitflag]=fmincon(@objective,y0,[],[],[],[],lower(active)/10000,upper(active)/10000,@constraints,op);
 end
catch err
 out.failure_reason=err.message;
end
out.evaluations=evals;out.elapsed_s=toc(start);
 function pl=decode(y)
  v=base;v(active)=y*10000;pl=plan;pl.times=v(1:N);pl.waits=v(N+1:end);
 end
 function evaluate(y)
  if isequal(cacheY,y),return;end
  evals=evals+1;cacheY=y;pl=decode(y);[cacheS,cacheR]=ctocscreen.rebuildArcPlan(pl,p,cfg);
  if cacheR.passed&&cacheR.total_dv_km_s<bestCost
   bestCost=cacheR.total_dv_km_s;out.plan=pl;out.schedule=cacheS;out.evaluation=cacheR;
  end
 end
 function v=objective(y)
  evaluate(y);v=1e5;if cacheR.completed,v=cacheR.total_dv_km_s;end
 end
 function v=residual(y)
  evaluate(y);v=ones(3*numel(intermediate),1)*100;
  if cacheR.completed,v=reshape(cacheR.position_residual_km(intermediate,:)',[],1)/1000;end
 end
 function [c,ceq]=constraints(y)
  evaluate(y);pl=decode(y);prev=[0;pl.times(ends(1:end-1))];
  c=[(prev+pl.waits+30-pl.times(begins))/10000;(30-diff([0;pl.times]))/10000];
  h=1e3;ceq=ones(3*numel(intermediate),1)*100;
  if cacheR.completed,h=(200-cacheR.min_altitude_km)/1000;ceq=reshape(cacheR.position_residual_km(intermediate,:)',[],1)/1000;end
  c=[c;h];
  if isfield(cfg,'visit_radius_km')
   miss=100*ones(numel(intermediate),1);
   if cacheR.completed,miss=cacheR.event_distances_km(intermediate).^2/cfg.search_visit_radius_km^2-1;end
   c=[c;miss];ceq=[];
  end
 end
 function yes=stop(~,~,~),yes=toc(start)>=cfg.window_seconds;end
 function yes=stopFeas(~,~,~),yes=toc(start)>=cfg.window_seconds*cfg.feasibility_fraction;end
end
