function report=eval_v4_endgame(files)
%EVAL_V4_ENDGAME Which targets make the expensive legs, and does H see them?
% Reads saved construction/search trajectories (diagnostic input only).
% For each leg (a burn followed by the target whose witness ends that leg),
% records the actual burn magnitude, the target's orbit class and the
% crossing-event estimate the node had for that target just before the burn.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(); mu=eph.model.mu;
cls=strings(35,1); cls([1 2 3])="GEO-eq"; cls(4:11)="GEO-incl"; cls(12:28)="GNSS";
cls([29 30])="MEO10k"; cls(31)="MEO8k-eq"; cls(32:34)="HEO"; cls(35)="HEO-GEO";
ecc=zeros(35,1); for j=1:35, x=eph.states0(j,:); r=x(1:3); v=x(4:6);
 ecc(j)=norm(((dot(v,v)-mu/norm(r))*r-dot(r,v)*v)/mu); end
rows=[]; names={};
for f=1:numel(files)
 s=load(fullfile(sim,files{f})); q=getq(s);
 [a,tr]=ctocscreen.v4.replay(q,eph,c);
 vis=find(a.distance_km<=1); [wt,o]=sort(a.witness_times_s(vis)); vis=vis(o);
 for k=1:numel(q.tau)
  next=find(wt>q.tau(k),1); if isempty(next), continue; end
  target=vis(next);
  % Node state just before this burn: truncate controls and replay.
  p=q; p.T=q.tau(k); p.tau=q.tau(1:k-1); p.u=q.u(1:k-1,:); p.witness(p.witness>p.T)=NaN;
  [pa,pt]=ctocscreen.v4.replay(p,eph,c); node=struct('q',pt.q,'actual',pa,'trace',pt);
  est=NaN; H=NaN;
  if pa.distance_km(target)>1
   [~,info]=ctocscreen.v4.events(node,eph,c); est=info.target_estimate(target); H=info.H;
  end
  rows(end+1,:)=[f,k,target,norm(q.u(k,:)),est,q.tau(k)/86400,ecc(target),H]; %#ok<AGROW>
 end
 names{f}=files{f};
end
T=array2table(rows,'VariableNames',{'file','burn','target','dv','estimate','day','target_ecc','H'});
T.class=cls(T.target);
expensive=T(T.dv>1,:);
fprintf('LEGS total=%d expensive(>1 km/s)=%d\n',height(T),height(expensive));
g=groupsummary(T,'class',{'mean','max'},{'dv','estimate'});
disp(g(:,{'class','GroupCount','mean_dv','max_dv','mean_estimate'}));
fprintf('expensive legs by class:\n'); disp(groupcounts(expensive,'class'));
fprintf('expensive legs day>=9: %d of %d\n',sum(expensive.day>=9),height(expensive));
ok=isfinite(T.estimate);
fprintf('corr(estimate,dv)=%.3f; expensive legs median estimate=%.3f vs actual median %.3f\n', ...
 corr(T.estimate(ok),T.dv(ok)),median(expensive.estimate,'omitnan'),median(expensive.dv));
fprintf('underestimate ratio actual/estimate: all median %.2f, expensive median %.2f\n', ...
 median(T.dv(ok)./T.estimate(ok)),median(expensive.dv./expensive.estimate,'omitnan'));
report=struct('diagnostic_only',true,'files',{names},'legs',T);
save(fullfile(sim,'runs/v4/development/schedule_eval/endgame.mat'),'report');
end
function q=getq(s)
if isfield(s,'report')&&isfield(s.report,'q'), q=s.report.q;
elseif isfield(s,'report')&&isfield(s.report,'best'), q=s.report.best.q;
elseif isfield(s,'result'), q=s.result.best.q;
else, error('Unknown file layout.'); end
end
