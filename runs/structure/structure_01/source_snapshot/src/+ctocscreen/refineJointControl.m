function out=refineJointControl(plan,focus,p,cfg,group)
%REFINEJOINTCONTROL Nested A/B/C continuous ablation, fixed branch labels.
% A: local times. B: same times + q. C/D: all times + q.
start=tic;N=numel(plan.ids);M=numel(plan.counts);
if ~isfield(plan,'aim_offsets_km'),plan.aim_offsets_km=zeros(M,3);end
if isfield(plan,'use_greedy_branches'),plan=rmfield(plan,'use_greedy_branches');end
[s,r]=ctocscreen.rebuildArcPlan(plan,p,cfg);
assert(r.passed,'Initial plan must be feasible.');
plan.locked_branch_ids=r.branch_ids;
ends=cumsum(plan.counts);begins=[1;ends(1:end-1)+1];
arcs=max(1,focus-1):min(M,focus+cfg.lookahead_arcs);
multi=find(plan.counts>1);
arcs=unique([arcs(:);multi;max(1,multi-1);min(M,multi+1)]);
events=[];for k=arcs',events=[events begins(k):ends(k)];end %#ok<AGROW>
events=unique(events);
if any(group=='CD'),arcs=(1:M)';events=1:N;end
aimIndices=reshape(1:3*M,3,M)';
active=[6+events(:);6+N+arcs(:);6+N+M+reshape(aimIndices(arcs,:)',[],1)];
if group~='A',active=[(1:6)';active];end
offset=[p.re_km+600;zeros(5+N+4*M,1)];
scale=[10;.001;.001;pi;2*pi;2*pi;10000*ones(N+M,1);ones(3*M,1)];
base=[plan.initial_q(:);plan.times;plan.waits;reshape(plan.aim_offsets_km',[],1)];x0=(base-offset)./scale;
lo=[p.re_km+590;-.001;-.001;0;plan.initial_q(5:6)'-pi;ones(N,1);zeros(M,1);-ones(3*M,1)];
hi=[p.re_km+610;.001;.001;pi;plan.initial_q(5:6)'+pi;p.horizon_s*ones(N+M,1);ones(3*M,1)];
if any(group=='AB')
 gap=diff([0;plan.times;p.horizon_s+3600]);radius=cfg.time_radius*min(gap(1:N),gap(2:N+1));
 lo(7:6+N)=max(1,plan.times-radius);hi(7:6+N)=min(p.horizon_s,plan.times+radius);
end
lo=(lo-offset)./scale;hi=(hi-offset)./scale;
% Linear time constraints, in units of 10000 s, include waits and chronology.
L=zeros(N+M,6+N+4*M);rhs=-cfg.minimum_gap_s/10000*ones(N+M,1);
for j=1:N,L(j,6+j)=-1;if j>1,L(j,5+j)=1;end,end
for m=1:M
 L(N+m,6+N+m)=1;L(N+m,6+begins(m))=-1;
 if m>1,L(N+m,6+ends(m-1))=1;end
end
fixed=setdiff(1:numel(x0),active);intermediate=setdiff(1:N,ends);
out=struct('plan',plan,'schedule',s,'evaluation',r,'initial_cost',r.total_dv_km_s, ...
 'evaluations',0,'invalid_evaluations',0,'exitflag',NaN,'output',struct(), ...
 'group',group,'active_indices',active,'elapsed_s',0,'failure_reason','');
cache=[];rr=[];ss=[];best=r.total_dv_km_s;
opt=optimoptions('fmincon','Algorithm','sqp','Display','off','UseParallel',false, ...
 'MaxIterations',cfg.joint_iterations,'MaxFunctionEvaluations',cfg.joint_evaluations, ...
 'FiniteDifferenceStepSize',1e-6,'ConstraintTolerance',1e-9, ...
 'OptimalityTolerance',1e-6,'StepTolerance',1e-9,'OutputFcn',@stop);
try
 [~,~,out.exitflag,out.output]=fmincon(@objective,x0(active),L(:,active), ...
  rhs-L(:,fixed)*x0(fixed),[],[],lo(active),hi(active),@constraints,opt);
catch err
 out.failure_reason=err.message;
end
out.elapsed_s=toc(start);
 function pl=decode(y)
  x=x0;x(active)=y;v=x.*scale+offset;pl=plan;
  pl.initial_q=v(1:6)';pl.initial_q(5:6)=mod(pl.initial_q(5:6),2*pi);
  pl.times=v(7:6+N);pl.waits=v(7+N:6+N+M);
  pl.aim_offsets_km=reshape(v(7+N+M:end),3,M)';
 end
 function evaluate(y)
  if isequal(cache,y),return;end
  cache=y;out.evaluations=out.evaluations+1;pl=decode(y);
  [ss,rr]=ctocscreen.rebuildArcPlan(pl,p,cfg);
  if ~rr.passed,out.invalid_evaluations=out.invalid_evaluations+1;end
  if rr.passed&&max(rr.event_distances_km)<=cfg.search_visit_radius_km&&rr.total_dv_km_s<best
   best=rr.total_dv_km_s;out.plan=pl;out.schedule=ss;out.evaluation=rr;
  end
 end
 function f=objective(y)
  evaluate(y);f=1e5;if rr.completed,f=rr.total_dv_km_s;end
 end
 function [c,eq]=constraints(y)
  evaluate(y);pl=decode(y);c=[(hypot(pl.initial_q(2),pl.initial_q(3))/.001)^2-(1-1e-9)^2;1e3];
  eq=[];miss=100*ones(numel(intermediate),1);
  if rr.completed
   c(2)=(200-rr.min_altitude_km)/1000;
   miss=rr.event_distances_km(intermediate).^2/cfg.search_visit_radius_km^2-1;
  end
  c=[c;sum(pl.aim_offsets_km.^2,2)/cfg.search_visit_radius_km^2-1;miss];
 end
 function yes=stop(~,~,~)
  yes=toc(start)>=cfg.chunk_seconds||isfile(cfg.stop_file);
 end
end
