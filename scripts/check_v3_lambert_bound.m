function check_v3_lambert_bound()
%CHECK_V3LAMBERTBOUND Compare geometric screening against broad enumeration.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
mu=398600.4418; limit=6; checked=0; oldIntervals=0; newIntervals=0;
for radius=[7000 26560 42164]
 for angle=[.3 1.7 3.0]
  a=[radius 0 0]; b=1.1*radius*[cos(angle) sin(angle) 0];
  for dt=[7200 40000 86400]
   full=ctocscreen.enumerateBranches(a,b,dt,mu,struct('max_revolutions',limit));
   n=ctocscreen.v3LambertRevolutions(a,b,dt,mu,limit);
   screened=ctocscreen.enumerateBranches(a,b,dt,mu,struct('max_revolutions',n));
   assert(numel(full)==numel(screened),'Lost a Lambert branch.');
   if ~isempty(full)
    assert(max(abs([full.v_depart]-[screened.v_depart]),[],'all')<1e-10);
   end
   checked=checked+1; oldIntervals=oldIntervals+limit; newIntervals=newIntervals+n;
  end
 end
end
fprintf('LAMBERT_BOUND: %d geometries/times preserve branches; multi-rev intervals %d -> %d.\n', ...
 checked,oldIntervals,newIntervals);
end
