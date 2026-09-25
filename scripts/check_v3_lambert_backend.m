function report=check_v3_lambert_backend()
%CHECK_V3_LAMBERT_BACKEND Compare branches, endpoints and runtime independently.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
mu=398600.4418; policy=struct('max_revolutions',6,'endpoint_tol_km',.001);
report=struct('cases',0,'old_seconds',0,'new_seconds',0,'old_branches',0, ...
 'new_branches',0,'fallbacks',0,'max_velocity_error_km_s',0);
for radius=[7000 26560 42164]
 for angle=[.01 .3 1.7 3.0 pi-1e-5]
  a=[radius 0 0]; b=1.1*radius*[cos(angle) sin(angle) 0];
  for dt=[3600 7200 40000 86400]
   local=policy; local.max_revolutions=ctocscreen.v3LambertRevolutions(a,b,dt,mu,policy.max_revolutions);
   clock=tic; old=ctocscreen.enumerateBranches(a,b,dt,mu,local); report.old_seconds=report.old_seconds+toc(clock);
   clock=tic; [fast,stats]=ctocscreen.v3LambertBranches(a,b,dt,mu,local); report.new_seconds=report.new_seconds+toc(clock);
   for k=1:numel(old)
    errors=arrayfun(@(v)norm(v.v_depart-old(k).v_depart)+norm(v.v_arrive-old(k).v_arrive),fast);
    assert(~isempty(errors)&&min(errors)<1e-6,'ctocscreen:v3:lambertComparison', ...
     'Missing legacy branch at radius %.0f angle %.9g dt %.0f.',radius,angle,dt);
    report.max_velocity_error_km_s=max(report.max_velocity_error_km_s,min(errors));
   end
   report.cases=report.cases+1; report.old_branches=report.old_branches+numel(old);
   report.new_branches=report.new_branches+numel(fast); report.fallbacks=report.fallbacks+stats.fallback;
  end
 end
end
disp(report);
end
