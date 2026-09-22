function out=refineStructureWindow(plan,core,p,cfg)
%REFINESTRUCTUREWINDOW Free q + incoming/core/outgoing window, whole-tour cost.
% Zero-visit burns are explicit free impulses; visit arcs use Lambert anchors.
start=tic;N=numel(plan.ids);M=numel(plan.counts);ends=cumsum(plan.counts);
positive=find(plan.counts>0);zero=find(plan.counts==0);intermediate=setdiff(1:N,ends(positive));
arcs=unique([1 max(1,core.first_arc-1):min(M,core.last_arc+cfg.lookahead_arcs)]');
arcs=unique([arcs(:);find(plan.counts>1);zero]);
events=[];for a=arcs',events=[events ends(a)-plan.counts(a)+1:ends(a)];end %#ok<AGROW>
events=unique(events);aim=reshape(1:3*M,3,M)';
active=[(1:6)';6+events(:);6+N+arcs;6+N+M+reshape(aim(intersect(arcs,positive),:)',[],1); ...
 6+N+4*M+(1:3*numel(zero))'];
scale=[10;.001;.001;pi;2*pi;2*pi;10000*ones(N+M,1);ones(3*M+3*numel(zero),1)];
offset=[p.re_km+600;zeros(numel(scale)-1,1)];
base=[plan.initial_q(:);plan.times;plan.waits;reshape(plan.aim_offsets_km',[],1);reshape(plan.zero_delta_v(zero,:)',[],1)];
x0=(base-offset)./scale;
lo=[p.re_km+590;-.001;-.001;0;plan.initial_q(5:6)'-pi;ones(N,1);zeros(M,1);-ones(3*M,1);base(7+N+4*M:end)-cfg.zero_velocity_window];
hi=[p.re_km+610;.001;.001;pi;plan.initial_q(5:6)'+pi;p.horizon_s*ones(N+M,1);ones(3*M,1);base(7+N+4*M:end)+cfg.zero_velocity_window];
gap=diff([0;plan.times;p.horizon_s+3600]);radius=cfg.time_radius*min(gap(1:N),gap(2:N+1));
lo(7:6+N)=max(1,plan.times-radius);hi(7:6+N)=min(p.horizon_s,plan.times+radius);
lo=(lo-offset)./scale;hi=(hi-offset)./scale;
L=zeros(N+M,numel(base));rhs=-cfg.minimum_gap_s/10000*ones(N+M,1);
for j=1:N,L(j,6+j)=-1;if j>1,L(j,5+j)=1;end,end
for m=1:M
 prior=find(plan.counts(1:m-1)>0,1,'last');if isempty(prior),lastEvent=0;firstArc=1;else,lastEvent=ends(prior);firstArc=prior+1;end
 L(N+m,6+N+(firstArc:m))=1;L(N+m,6+lastEvent+1)=-1;
 if lastEvent>0,L(N+m,6+lastEvent)=1;end
 if m>1&&plan.counts(m-1)==0,lo(6+N+m)=cfg.minimum_gap_s/10000;end
end
fixed=setdiff(1:numel(base),active);cache=[];rr=[];ss=[];best=Inf;phase='feasibility';bestY=x0(active);
out=struct('plan',plan,'schedule',[],'evaluation',struct('passed',false), ...
 'evaluations',0,'invalid_evaluations',0,'exitflag',NaN,'feasibility_exitflag',NaN, ...
 'active_indices',active,'failure_reason','','elapsed_s',0);
opt=optimoptions('fmincon','Algorithm','sqp','Display','off','UseParallel',false, ...
 'MaxIterations',150,'MaxFunctionEvaluations',cfg.joint_evaluations,'FiniteDifferenceStepSize',1e-6, ...
 'ConstraintTolerance',1e-8,'OptimalityTolerance',1e-6,'StepTolerance',1e-9,'OutputFcn',@stop);
try
 evaluate(x0(active));
 if ~rr.completed,out.failure_reason='Initial branch construction failed';out.elapsed_s=toc(start);return;end
 if ~out.evaluation.passed
  [y,~,out.feasibility_exitflag,out.feasibility_output]=fmincon(@objective,x0(active),L(:,active), ...
   rhs-L(:,fixed)*x0(fixed),[],[],lo(active),hi(active),@constraints,opt);
  if ~out.evaluation.passed,bestY=y;end
 end
 phase='cost';cache=[];
 if toc(start)<cfg.chunk_seconds
  [~,~,out.exitflag,out.output]=fmincon(@objective,bestY,L(:,active),rhs-L(:,fixed)*x0(fixed), ...
   [],[],lo(active),hi(active),@constraints,opt);
 end
catch err,out.failure_reason=err.message;end
out.elapsed_s=toc(start);
 function pl=decode(y)
  xx=x0;xx(active)=y;v=xx.*scale+offset;pl=plan;
  pl.initial_q=v(1:6)';pl.initial_q(5:6)=mod(pl.initial_q(5:6),2*pi);
  pl.times=v(7:6+N);pl.waits=v(7+N:6+N+M);
  pl.aim_offsets_km=reshape(v(7+N+M:6+N+4*M),3,M)';
  pl.zero_delta_v(zero,:)=reshape(v(7+N+4*M:end),3,[])';
 end
 function evaluate(y)
  if isequal(cache,y),return;end
  cache=y;out.evaluations=out.evaluations+1;pl=decode(y);[ss,rr]=ctocscreen.rebuildArcPlan(pl,p,cfg);
  if ~rr.passed,out.invalid_evaluations=out.invalid_evaluations+1;end
  if rr.passed&&max(rr.event_distances_km)<=cfg.search_visit_radius_km+1e-6&&rr.total_dv_km_s<best
   best=rr.total_dv_km_s;bestY=y;out.plan=pl;out.schedule=ss;out.evaluation=rr;
   u=vecnorm(ss.delta_v_km_s,2,2);
   out.cost_parts_km_s=[sum(u(1:core.first_arc-1)) sum(u(core.first_arc:core.last_arc)) sum(u(core.last_arc+1:end))];
  end
 end
 function f=objective(y)
  evaluate(y);f=1e5;
  if rr.completed
   if strcmp(phase,'feasibility'),f=sum(max(0,rr.event_distances_km-cfg.search_visit_radius_km).^2)/1e6;
   else,f=rr.total_dv_km_s;end
  end
 end
 function [c,eq]=constraints(y)
  evaluate(y);pl=decode(y);eq=[];
  c=[(hypot(pl.initial_q(2),pl.initial_q(3))/.001)^2-(1-1e-9)^2; ...
   sum(pl.aim_offsets_km(positive,:).^2,2)/cfg.search_visit_radius_km^2-1;1e3];
  if rr.completed,c(end)=(200-rr.min_altitude_km)/1000;end
  if strcmp(phase,'cost')
   miss=100*ones(numel(intermediate),1);
   if rr.completed,miss=(rr.event_distances_km(intermediate)-cfg.search_visit_radius_km)/1000;end
   c=[c;miss];
  end
 end
 function yes=stop(~,~,~)
  limit=cfg.chunk_seconds;if strcmp(phase,'feasibility'),limit=cfg.chunk_seconds*.6;end
  yes=toc(start)>=limit||isfile(cfg.stop_file);
 end
end
