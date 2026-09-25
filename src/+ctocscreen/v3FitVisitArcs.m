function [s,report]=v3FitVisitArcs(s,eph,c,budget)
%V3FITVISITARCS Physical initial guess fitting every planned visit on each arc.
% Overdetermined arcs may remain infeasible; the complete joint solve follows.
clock=tic; m=eph.model; M=numel(s.maneuver_times_s);
report=struct('passed',false,'reason','','elapsed_s',0,'arcs',{{}}, ...
 'method','all_visit_arc_least_squares','planned_target_ids',find(isfinite(s.visit_plan_times_s)).');
fine=c; fine.shooting_reltol=min(c.shooting_reltol,1e-12);
try
 x=ctocscreen.initialState(s.initial_q,m.mu,m.re).';
 if M>0, first=s.maneuver_times_s(1); else, first=s.duration_s; end
 x=ctocscreen.v3Arc(x,0,first,m,fine);
 ends=[s.maneuver_times_s(2:end);s.duration_s];
 for k=1:M
  assert(toc(clock)<budget,'ctocscreen:v3:visitFitBudget','Visit construction budget exhausted.');
  ta=s.maneuver_times_s(k); tb=ends(k); v=x(4:6)+s.delta_v_km_s(k,:).';
  ids=find(s.visit_plan_times_s>ta & s.visit_plan_times_s<=tb);
  times=s.visit_plan_times_s(ids);
  item=struct('target_ids',ids.','before_km',[],'after_km',[],'iterations',0,'branch_reseeded',false);
  if ~isempty(ids)
   goals=ctocscreen.v3QueryTargets(eph,ids,times);
   try
    [res,B]=residual(v);
   catch
    [last,lastId]=max(times); remaining=c; remaining.budget_s=max(0,budget-toc(clock));
    dv=ctocscreen.v3GuidedTransfer(x,ta,last,goals(lastId,:).',m,remaining,tic);
    v=x(4:6)+dv; item.branch_reseeded=true; [res,B]=residual(v);
   end
   item.before_km=vecnorm(reshape(res,3,[]),2,1);
   for it=1:c.connection_iterations
    if toc(clock)>=budget||max(vecnorm(reshape(res,3,[]),2,1))<=c.connection_tolerance_km, break; end
    step=[B;1e-5*eye(3)]\[-res;zeros(3,1)];
    step=step*min(1,2/max(norm(step),eps)); accepted=false;
    for alpha=[1 .5 .25 .125]
     if toc(clock)>=budget, break; end
     try
      [rr,BB]=residual(v+alpha*step);
      if norm(rr)<norm(res), v=v+alpha*step; res=rr; B=BB; accepted=true; break; end
     catch
      % Failed propagation is an invalid proposal, not a finite penalty.
     end
    end
    item.iterations=it;
    if ~accepted, break; end
   end
   item.after_km=vecnorm(reshape(res,3,[]),2,1);
  end
  s.delta_v_km_s(k,:)=(v-x(4:6)).';
  x=ctocscreen.v3Arc([x(1:3);v],ta,tb,m,fine);
  report.arcs{end+1}=item;
 end
 s=ctocscreen.v3Normalize(s,m); report.passed=true;
catch err
 report.reason=[err.identifier ': ' err.message];
end
report.elapsed_s=toc(clock);
 function [r,B]=residual(velocity)
  [~,~,sol]=ctocscreen.v3Arc([x(1:3);velocity],ta,max(times),m,fine,true);
  y=deval(sol,times); r=reshape(y(1:3,:)-goals.',[],1); B=zeros(3*numel(ids),3);
  for sample=1:numel(ids)
   P=reshape(y(7:end,sample),6,6); B(3*(sample-1)+(1:3),:)=P(1:3,4:6);
  end
 end
end
