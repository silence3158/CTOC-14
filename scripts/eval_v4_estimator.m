function report=eval_v4_estimator(files)
%EVAL_V4_ESTIMATOR Does a Lambert-based leg estimate track realized legs?
% For each leg of saved chains (diagnostic input only), compares the realized
% burn magnitude with (a) events.m two-body tangential estimate and (b) the
% cheapest crossing.m Lambert proposal to that target from the pre-burn node.
% Hypothesis: (b) correlates better, so it should replace (a) in H.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(struct('event_candidates',200,'events_per_target',3)); stream=RandStream('mt19937ar','Seed',1);
rows=zeros(0,6);
for f=1:numel(files)
 s=load(fullfile(sim,files{f})); q=s.report.q;
 [a,~]=ctocscreen.v4.replay(q,eph,c);
 vis=find(a.distance_km<=1); [wt,o]=sort(a.witness_times_s(vis)); vis=vis(o);
 for k=1:numel(q.tau)
  next=find(wt>q.tau(k),1); if isempty(next), continue; end
  target=vis(next);
  % Node at the previous arrival (the burn may be delayed along its coast).
  prev=find(wt<=q.tau(k),1,'last'); tEnd=0; if ~isempty(prev), tEnd=wt(prev); end
  p=q; keep=p.tau<tEnd|(p.tau<=tEnd&false); p.tau=q.tau(keep); p.u=q.u(keep,:); p.T=max(tEnd,max([p.tau;0]));
  p.witness(p.witness>p.T)=NaN;
  [pa,pt]=ctocscreen.v4.replay(p,eph,c); if pa.distance_km(target)<=1, continue; end
  node=struct('q',pt.q,'actual',pa,'trace',pt);
  [ev,info]=ctocscreen.v4.events(node,eph,c);
  e1=info.target_estimate(target);
  mask=[ev.id]==target; lam=NaN;
  if any(mask)
   cand=ctocscreen.v4.crossing(node,eph,c,ev(mask),info,stream,20);
   if ~isempty(cand), lam=min([cand.dv]); end
  end
  rows(end+1,:)=[f,k,target,norm(q.u(k,:)),e1,lam]; %#ok<AGROW>
 end
end
ok=all(isfinite(rows(:,4:6)),2); R=rows(ok,:);
c1=corr(R(:,5),R(:,4)); c2=corr(R(:,6),R(:,4));
s1=corr(R(:,5),R(:,4),'type','Spearman'); s2=corr(R(:,6),R(:,4),'type','Spearman');
fprintf('ESTIMATOR legs=%d usable=%d\n',size(rows,1),size(R,1));
fprintf('two-body tangential: Pearson %.3f Spearman %.3f median actual/est %.2f\n',c1,s1,median(R(:,4)./R(:,5)));
fprintf('Lambert crossing   : Pearson %.3f Spearman %.3f median actual/est %.2f\n',c2,s2,median(R(:,4)./R(:,6)));
big=R(:,4)>1; fprintf('legs >1 km/s: n=%d, median est tangential %.3f, Lambert %.3f, actual %.3f\n', ...
 sum(big),median(R(big,5)),median(R(big,6)),median(R(big,4)));
report=struct('diagnostic_only',true,'rows',rows,'pearson',[c1 c2],'spearman',[s1 s2]);
save(fullfile(sim,'runs/v4/development/schedule_eval/estimator.mat'),'report');
end
