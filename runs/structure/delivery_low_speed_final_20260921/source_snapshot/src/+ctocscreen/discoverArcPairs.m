function candidates=discoverArcPairs(s,r,p,cfg,stream)
%DISCOVERARCPAIRS Rank all future targets by continuation geometry, not ID proximity.
% Coarse discovery is only a proposal heuristic, never an encounter certificate.
N=numel(s.event_times_s);candidates=zeros(0,5);
anchors=randperm(stream,N-1);anchors=anchors(1:min(cfg.discovery_anchors,N-1));
for a=anchors
 lo=s.event_times_s(a)+60;
 hi=min(p.horizon_s,s.event_times_s(min(N,a+2))-60);if hi<=lo,continue;end
 grid=linspace(lo,hi,cfg.discovery_samples);best=inf(N-a,1);when=best;
 for tt=grid
  x=ctocscreen.propagateTwoBody(r.event_states(a,:),tt-s.event_times_s(a),p.mu_km3_s2);
  targets=ctocscreen.targetStates(p,s.event_target_ids(a+1:end),tt);
  distance=vecnorm(targets(:,1:3)-x(1:3),2,2);improve=distance<best;
  best(improve)=distance(improve);when(improve)=tt;
 end
 for b=(a+1):N
  dt=when(b-a)-s.event_times_s(a);
  % Velocity-scale correction proxy, plus relocation burden. No physical constraint.
  score=best(b-a)/max(600,dt)+cfg.relocation_weight*(b-a-1);
  candidates(end+1,:)=[a b when(b-a) best(b-a) score]; %#ok<AGROW>
 end
end
[~,order]=sort(candidates(:,5));candidates=candidates(order,:);
end
