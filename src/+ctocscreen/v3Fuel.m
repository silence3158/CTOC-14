function [best,raw,report]=v3Fuel(p,z,c,budget)
%V3FUEL SQP local steps, preserving best feasible iterate and raw diagnostics.
if strcmp(c.fuel_algorithm,'sparse_sqp')
 [best,raw,report]=ctocscreen.v3SparseFuel(p,z,c,budget); return
end
clock=tic; best=z; raw=z; start=z; bestCost=p.objective(z); initial=[]; failure='';
opt=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
 'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, ...
 'ConstraintTolerance',min(c.constraint_tolerance,1e-11), ...
 'MaxIterations',c.joint_iterations,'MaxFunctionEvaluations',max(100,10*c.joint_iterations), ...
 'OutputFcn',@observe);
flag=0; output=struct(); iterations=0;
try
 [raw,~,flag,output]=fmincon(p.objective,z,full(p.A),p.b,[],[], ...
  max(p.lb,z-c.fuel_trust_radius),min(p.ub,z+c.fuel_trust_radius),@denseConstraints,opt);
 raw=min(p.ub,max(p.lb,raw)); consider(raw);
catch err
 failure=[err.identifier ': ' err.message];
end
[ci,eq]=p.constraints(raw);
report=struct('algorithm','sqp','exitflag',flag,'output',output,'initial',initial, ...
 'raw_violation',max([0;ci;abs(eq);p.A*raw-p.b;p.lb-raw;raw-p.ub]), ...
 'raw_objective',p.objective(raw),'selected_objective',bestCost, ...
 'failure',failure,'elapsed_s',toc(clock),'trust_radius',c.fuel_trust_radius);
report.refinement_performed=iterations>0||flag>0||flag==-2;
 function [ci,eq,G,GE]=denseConstraints(x)
  [ci,eq,G,GE]=p.constraints(x); G=full(G); GE=full(GE);
 end
 function consider(x)
  [ci,eq]=p.constraints(x);
  if max([0;ci;abs(eq);p.A*x-p.b;p.lb-x;x-p.ub])<=min(10*c.constraint_tolerance,1e-10) && p.objective(x)<bestCost
   best=x; bestCost=p.objective(x);
  end
 end
 function stop=observe(x,values,phase)
  raw=x; consider(x);
  iterations=max(iterations,values.iteration);
  if strcmp(phase,'init')
   [ci,eq]=p.constraints(x); initial=struct('change',norm(x-start),'violation',max([0;ci;abs(eq)]));
  end
  stop=toc(clock)>=budget||(~isempty(c.stop_file)&&isfile(c.stop_file));
 end
end
