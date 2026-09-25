function [z,report]=v3Restore(p,c,budget,start)
%V3RESTORE Damped Gauss-Newton with hard linear/bound constraints and backtracking.
% All joint variables participate. Merit is NOT a trajectory delta-V.
if nargin<4, start=p.z0; end
clock=tic; z=start; n=numel(z); trust=.1; damping=1e-6; rejected=0;
[ci,eq,Gi,Ge]=p.constraints(z); residual=[max(ci,0);eq]; merit=.5*sum(residual.^2);
initial=merit; history=merit; flags=[]; calls=1; it=0; evaluatedSteps=0; solverWork=false;
opt=optimoptions('lsqlin','Display','off','Algorithm','interior-point', ...
 'MaxIterations',100,'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);
while it<c.restoration_iterations && toc(clock)<budget
 if max([0;ci;abs(eq);p.A*z-p.b])<=c.constraint_tolerance, break; end
 if ~isempty(c.stop_file)&&isfile(c.stop_file), break; end
 it=it+1; active=ci>0;
 J=[Gi(:,active).';Ge.']; r=[ci(active);eq];
 matrix=[J;sqrt(damping)*speye(n)]; right=[-r;zeros(n,1)];
 lo=max(p.lb-z,-trust); hi=min(p.ub-z,trust);
 [step,~,~,flag,output]=lsqlin(matrix,right,p.A,p.b-p.A*z,[],[],lo,hi,[],opt);
 solverWork=solverWork||output.iterations>0||flag>0||flag==-2;
 flags(end+1)=flag; %#ok<AGROW>
 if isempty(step)||any(~isfinite(step)), damping=damping*10; trust=trust/2; continue; end
 accepted=false;
 for alpha=[1 .5 .25 .125 .0625 .03125]
  if toc(clock)>=budget, break; end
  trial=z+alpha*step;
  % Preserve the formal eccentricity disk during restoration. Never decode
  % an out-of-disk iterate and call the ensuing rejection a solver result.
  q=p.q_indices(2:3); radius=norm(trial(q));
  if radius>1-1e-8, trial(q)=trial(q)*(1-1e-8)/radius; end
  try
   [tc,te,tG,tGe]=p.constraints(trial); calls=calls+1;
   evaluatedSteps=evaluatedSteps+all(isfinite([tc;te]));
   f=.5*(sum(max(tc,0).^2)+sum(te.^2));
   if isfinite(f)&&f<merit
    z=trial; ci=tc; eq=te; Gi=tG; Ge=tGe; merit=f; accepted=true; break
   end
  catch
   rejected=rejected+1; % Invalid integration is rejected, never a finite cost.
  end
 end
 history(end+1)=merit; %#ok<AGROW>
 if accepted, trust=min(1,trust*1.5); damping=max(1e-12,damping/3);
 else, trust=trust/2; damping=damping*10; end
 if trust<1e-9, break; end
end
violation=max([0;ci;abs(eq);p.A*z-p.b]); exitflag=0;
if violation<=c.constraint_tolerance, exitflag=1; end
reason='stop_requested';
if exitflag==1, reason='constraints_restored';
elseif it>=c.restoration_iterations, reason='iteration_limit';
elseif toc(clock)>=budget, reason='time_budget';
elseif trust<1e-9, reason='trust_region_stalled'; end
report=struct('initial_merit',initial,'final_merit',merit,'merit_history',history, ...
 'exitflag',exitflag,'output',struct('iterations',it,'funcCount',calls,'subproblem_exitflags',flags), ...
 'elapsed_s',toc(clock),'violation',violation,'invalid_trials',rejected, ...
 'refinement_performed',solverWork||evaluatedSteps>0,'evaluated_steps',evaluatedSteps, ...
 'method','damped joint Gauss-Newton / constrained lsqlin / backtracking','termination_reason',reason);
end
