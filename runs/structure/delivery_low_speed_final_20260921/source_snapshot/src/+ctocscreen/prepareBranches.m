function [c,dp]=prepareBranches(c,p,cfg)
%PREPAREBRANCHES Nominal branch DP followed separately by actual-state replay.
s=ctocscreen.initialState(c.initial_q,p.mu_km3_s2,p.re_km);
[s,~]=ctocscreen.propagateTwoBody(s,c.wait_s,p.mu_km3_s2);
v0=s(4:6); t=c.wait_s; sets=cell(1,35);
for k=1:35
 t=t+c.tof_s(k); [target,~]=ctocscreen.targetStates(p,c.order(k),t);
 [b,~]=ctocscreen.enumerateBranches(s(1:3),target(1:3),c.tof_s(k),p.mu_km3_s2,cfg.branch_policy);
 keep=false(size(b));
 for j=1:numel(b)
  a=ctocscreen.checkArc([s(1:3) b(j).v_depart],c.tof_s(k),p);
  keep(j)=strcmp(a.status,'ok') && a.min_altitude_km>=200;
 end
 sets{k}=b(keep); s=target;
end
dp=ctocscreen.chooseBranchesDP(sets,v0);
if strcmp(dp.status,'ok'), c.branch_ids=dp.branch_ids; c.branch_refs=dp.selected; end
end
