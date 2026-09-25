function [raw,report]=v3JointDirection(p,radius)
%V3JOINTDIRECTION One sparse joint QP proposal, not a feasibility claim.
[ci,eq,G,GE]=p.constraints(p.z0); [f,g]=p.objective(p.z0);
[Av,bv]=p.encounter_linearization(p.z0);
opt=optimoptions('quadprog','Display','off','Algorithm','interior-point-convex', ...
 'MaxIterations',60,'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-9);
[step,~,flag,output]=quadprog(speye(numel(p.z0)),g,[p.A;G.';Av], ...
 [p.b-p.A*p.z0;-ci;bv],GE.',-eq,max(p.lb-p.z0,-radius),min(p.ub-p.z0,radius),[],opt);
raw=p.z0;
if flag>0 && all(isfinite(step)), raw=raw+step; end
report=struct('initial',struct('change',0),'initial_objective',f, ...
 'raw_objective',p.objective(raw),'exitflag',flag,'output',output,'radius',radius, ...
 'refinement_performed',output.iterations>0||flag>0||flag==-2);
end
