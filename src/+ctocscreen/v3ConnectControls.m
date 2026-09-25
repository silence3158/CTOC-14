function [s,report]=v3ConnectControls(s,anchors,m,c,budget)
%V3CONNECTCONTROLS Retract a proposal onto continuous, real-state arcs.
% anchors(k,:) is the requested position at the NEXT burn (or final time).
% Only existing impulses are corrected. Never reset a propagated state.
clock=tic; M=numel(s.maneuver_times_s);
report=struct('passed',false,'reason','','elapsed_s',0,'endpoint_errors_km',nan(M,1), ...
 'correction_km_s',zeros(M,1),'iterations',zeros(M,1),'branch_reseeds',0,'reused_integrations',0);
try
 assert(isequal(size(anchors),[M 3]),'ctocscreen:v3:anchors','Wrong anchor dimensions.');
 assert(all(diff(s.maneuver_times_s)>0),'ctocscreen:v3:anchors','Projection needs distinct ordered burns.');
 x=ctocscreen.initialState(s.initial_q,m.mu,m.re).';
 if M==0, report.passed=true; report.elapsed_s=toc(clock); return; end
 x=ctocscreen.v3Arc(x,0,s.maneuver_times_s(1),m,c,false,true);
 fine=c; fine.shooting_reltol=min(c.shooting_reltol,1e-13); fine.max_step_s=min(60,c.max_step_s);
 ends=[s.maneuver_times_s(2:end);s.duration_s];
 for k=1:M
  assert(toc(clock)<budget,'ctocscreen:v3:connectBudget','Connection budget exhausted.');
  ta=s.maneuver_times_s(k); tb=ends(k); v=x(4:6)+s.delta_v_km_s(k,:).'; original=v; reseeded=false; available=false;
  if tb==ta
   report.endpoint_errors_km(k)=norm(x(1:3)-anchors(k,:).'); continue
  end
  for it=1:c.connection_iterations
   assert(toc(clock)<budget,'ctocscreen:v3:connectBudget','Connection budget exhausted.');
   start=[x(1:3);v];
   if ~available
    [y,~,sol]=ctocscreen.v3Arc(start,ta,tb,m,fine,false,true); available=true;
   else
    report.reused_integrations=report.reused_integrations+1;
   end
   residual=y(1:3)-anchors(k,:).'; err=norm(residual); report.iterations(k)=it;
   if err<=c.connection_tolerance_km, break; end
   [~,P]=ctocscreen.v3Arc(start,ta,tb,m,fine,true);
   B=P(1:3,4:6);
   if rcond(B)<=1e-12 && ~reseeded
    v=reseed(); reseeded=true; available=false; continue
   end
   assert(rcond(B)>1e-12,'ctocscreen:v3:connectSingular','Singular arc position sensitivity.');
   step=B\residual; accepted=false;
   for alpha=[1 .5 .25 .125 .0625]
    [trial,~,trialSol]=ctocscreen.v3Arc([x(1:3);v-alpha*step],ta,tb,m,fine,false,true);
    if norm(trial(1:3)-anchors(k,:).')<err
     v=v-alpha*step; y=trial; sol=trialSol; available=true; accepted=true; break
    end
   end
   if ~accepted && ~reseeded
    v=reseed(); reseeded=true; available=false; continue
   end
   assert(accepted,'ctocscreen:v3:connectStalled','Arc correction did not reduce position error.');
  end
  if ~available
   [y,~,sol]=ctocscreen.v3Arc([x(1:3);v],ta,tb,m,fine,false,true);
  else
   report.reused_integrations=report.reused_integrations+1;
  end
  report.endpoint_errors_km(k)=norm(y(1:3)-anchors(k,:).');
  assert(report.endpoint_errors_km(k)<=c.connection_tolerance_km, ...
   'ctocscreen:v3:connectMiss','Arc correction did not converge.');
  height=ctocscreen.v3Height(sol,m,fine);
  assert(height.passed,'ctocscreen:v3:connectAltitude','Connected arc violates altitude gate.');
  s.delta_v_km_s(k,:)=(v-x(4:6)).'; report.correction_km_s(k)=norm(v-original);
  x=y; % Actual arrival position AND velocity, not the anchor.
 end
 s=ctocscreen.v3Normalize(s,m); report.passed=true;
catch err
 report.reason=[err.identifier ': ' err.message];
end
report.elapsed_s=toc(clock);
 function velocity=reseed()
  remaining=fine; remaining.budget_s=max(0,budget-toc(clock));
  [dv,~]=ctocscreen.v3GuidedTransfer(x,ta,tb,anchors(k,:).',m,remaining,tic);
  velocity=x(4:6)+dv; report.branch_reseeds=report.branch_reseeds+1;
 end
end
