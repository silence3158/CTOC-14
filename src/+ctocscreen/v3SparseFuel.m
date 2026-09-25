function [best,raw,report]=v3SparseFuel(p,z,c,budget)
%V3SPARSEFUEL Sparse regularized SQP steps without a dense BFGS Hessian.
clock=tic; start=z; best=z; raw=z; bestCost=p.objective(z); failure=''; flags=[]; evaluatedSteps=0; solverWork=false;
[ci,eq,G,GE]=p.constraints(z); [f,g]=p.objective(z);
initial=struct('change',0,'violation',max([0;ci;abs(eq)]));
opt=optimoptions('quadprog','Display','off','Algorithm','interior-point-convex', ...
 'MaxIterations',30,'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-9);
H=speye(numel(z)); radius=c.fuel_trust_radius;
for it=1:c.joint_iterations
 if toc(clock)>=budget||(~isempty(c.stop_file)&&isfile(c.stop_file)), break; end
 try
  [step,~,flag,output]=quadprog(H,g,[p.A;G.'],[p.b-p.A*z;-ci],GE.',-eq, ...
   max(p.lb-z,-radius),min(p.ub-z,radius),[],opt);
  solverWork=solverWork||output.iterations>0||flag>0||flag==-2;
  flags(end+1)=flag;
  if isempty(step)||any(~isfinite(step))||flag<=0, radius=radius/2; continue; end
  merit=f+1e4*(norm(eq,1)+sum(max(ci,0))); accepted=false;
  for alpha=[1 .3 .1 .03 .01]
   if toc(clock)>=budget, break; end
   trial=min(p.ub,max(p.lb,z+alpha*step));
   [tc,te,tG,tGE]=p.constraints(trial); [tf,tg]=p.objective(trial);
   evaluatedSteps=evaluatedSteps+all(isfinite([tc;te;tf]));
   if tf+1e4*(norm(te,1)+sum(max(tc,0)))<merit
    z=trial; raw=z; f=tf; g=tg; ci=tc; eq=te; G=tG; GE=tGE; accepted=true;
    if max([0;ci;abs(eq);p.A*z-p.b])<=min(10*c.constraint_tolerance,1e-10)&&f<bestCost
     best=z; bestCost=f;
    end
    break
   end
  end
  if ~accepted, radius=radius/2; end
  if radius<1e-9, break; end
 catch err
  failure=[err.identifier ': ' err.message]; break
 end
end
[ci,eq]=p.constraints(raw);
report=struct('algorithm','sparse_sqp','exitflag',0,'output',struct('iterations',numel(flags),'qp_flags',flags), ...
 'initial',initial,'raw_violation',max([0;ci;abs(eq);p.A*raw-p.b]), ...
 'raw_objective',p.objective(raw),'selected_objective',bestCost,'failure',failure, ...
 'refinement_performed',solverWork||evaluatedSteps>0,'evaluated_steps',evaluatedSteps, ...
 'elapsed_s',toc(clock),'trust_radius',c.fuel_trust_radius,'initial_objective',p.objective(start));
end
