function diagnose_v3_fuel()
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
a=load(fullfile(sim,'runs/v3/search/v3_budget_480_01/checkpoint.mat')); c=ctocscreen.v3Defaults(a.state.config);
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
entry=a.state.recovery_pool{1};
if isfield(entry,'schedule')
 s=entry.schedule; c.joint_target_ids=entry.target_ids;
else
 s=entry; c.joint_target_ids=find(isfinite(s.witness_times_s)).';
end
c.shooting_reltol=1e-12;
p=ctocscreen.v3ShootingProblem(s,eph,c); results={};
for algorithm={'sqp'}
 clock=tic; initial=[]; best=p.z0;
 opt=optimoptions('fmincon','Display','off','Algorithm',algorithm{1}, ...
  'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, ...
  'MaxIterations',60,'MaxFunctionEvaluations',1000,'ConstraintTolerance',1e-11,'OutputFcn',@stop);
 if strcmp(algorithm{1},'interior-point'), opt=optimoptions(opt,'HessianApproximation','lbfgs'); end
 [z,f,flag,out]=fmincon(p.objective,p.z0,full(p.A),p.b,[],[],max(p.lb,p.z0-1e-4),min(p.ub,p.z0+1e-4),p.constraints,opt);
 z=min(p.ub,max(p.lb,z));
 [ci,eq]=p.constraints(z); if max([0;ci;abs(eq)])>1e-10||p.objective(best)<p.objective(z), z=best; end; f=p.objective(z);
 [ci,eq]=p.constraints(z); r=ctocscreen.v3Replay(p.decode(z),eph,c,true);
 result=struct('algorithm',algorithm{1},'initial',initial,'flag',flag,'output',out,'f',f, ...
  'violation',max([0;ci;abs(eq)]),'replay',r,'elapsed_s',toc(clock)); results{end+1}=result;
 fprintf('%s cost %.9f -> %.9f flag %d violation %.9g visits %d height %d elapsed %.2f\n',algorithm{1},p.objective(p.z0),f,flag,result.violation,r.visit_count,r.height_passed,result.elapsed_s); disp(initial);
end
save(fullfile(sim,'runs/v3/development/fuel_diagnosis.mat'),'results');
 function yes=stop(z,~,phase)
  [ci,eq]=p.constraints(z);
  if max([0;ci;abs(eq);p.A*z-p.b;p.lb-z;z-p.ub])<=1e-10 && p.objective(z)<p.objective(best), best=z; end
  if strcmp(phase,'init')
   [ci,eq]=p.constraints(z); initial=struct('change',norm(z-p.z0),'violation',max([0;ci;abs(eq)]));
  end
  yes=toc(clock)>20;
 end
end
