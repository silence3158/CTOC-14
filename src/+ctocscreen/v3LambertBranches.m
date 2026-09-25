function [branches,report]=v3LambertBranches(r1,r2,dt,mu,policy)
%V3LAMBERTBRANCHES Vendored Izzo/Gooding seeds with endpoint checks and fallback.
report=struct('fallback',false,'checked',0,'extra_failures',0);
branches=struct('v_depart',{},'v_arrive',{},'rev_estimate',{},'geometry',{});
r1=r1(:).'; r2=r2(:).'; fallback=false;
count=ctocscreen.v3LambertRevolutions(r1,r2,dt,mu,policy.max_revolutions);
for direction=[1 -1]
 for revolutions=[0 reshape([1:count;-(1:count)],1,[])]
  try
   [v,w,~,flag]=lambert_izzo_gooding(r1,r2,direction*dt/86400,revolutions,mu);
   if flag==-1, continue; end
   if flag~=1||~isreal(v)||~isreal(w)||any(~isfinite([v w]))
    fallback=true; continue;
   end
   [endState,check]=ctocscreen.propagateTwoBody([r1 v],dt,mu);
   report.checked=report.checked+1;
   if ~strcmp(check.status,'ok')||norm(endState(1:3)-r2)>policy.endpoint_tol_km
    fallback=true; continue;
   end
   if any(arrayfun(@(b)norm(v-b.v_depart)<1e-9,branches)), continue; end
   geometry='short'; if direction<0, geometry='long'; end
   branches(end+1)=struct('v_depart',v,'v_arrive',w,'rev_estimate',abs(revolutions),'geometry',geometry);
  catch
   fallback=true; report.extra_failures=report.extra_failures+1;
  end
 end
end
if fallback
 legacy=ctocscreen.enumerateBranches(r1,r2,dt,mu,policy); report.fallback=true;
 for k=1:numel(legacy)
  b=legacy(k);
  if any(arrayfun(@(a)norm(a.v_depart-b.v_depart)<1e-9,branches)), continue; end
  branches(end+1)=struct('v_depart',b.v_depart,'v_arrive',b.v_arrive, ...
   'rev_estimate',b.rev_estimate,'geometry',char(b.geometry));
 end
end
end
