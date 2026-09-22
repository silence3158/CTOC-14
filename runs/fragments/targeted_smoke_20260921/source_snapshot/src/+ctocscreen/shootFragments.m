function [f,proposal,diagnostic]=shootFragments(plan,ev,p,cfg,stream)
%SHOOTFRAGMENTS Target-guided shot then discover the second target on the arc.
f=struct('passed',false); proposal=plan;
diagnostic=struct('shot',1,'gated',0,'refined',0,'mode','guided','k',0,'ids',[], ...
 'screen_distance_km',Inf);
n=numel(plan.ids); available=find(~plan.skip(1:n-1)&~plan.skip(2:n));
if isempty(available), return; end
k=available(randi(stream,numel(available))); diagnostic.k=k;
tin=0; xin=ctocscreen.initialState(plan.initial_q,p.mu_km3_s2);
if k>1, tin=plan.times(k-1); xin=ev.event_states(k-1,:); end
finish=p.horizon_s; if k+2<=n, finish=plan.times(k+2)-60; end
span=finish-tin; oldA=plan.times(k)-tin;
w=rand(stream)*min(.4*oldA,span*.2);
da=min(span*.7,max(300,oldA*(.7+.6*rand(stream))-w));
xp=ctocscreen.propagateTwoBody(xin,w,p.mu_km3_s2);
ta=tin+w+da; A=plan.ids(k); tar=ctocscreen.targetStates(p,A,ta);
[bs,~]=ctocscreen.enumerateBranches(xp(1:3),tar(1:3),da,p.mu_km3_s2,cfg.branch_policy);
if isempty(bs)||finish-ta<300, return; end
[~,ib]=min(arrayfun(@(b)norm(b.v_depart-xp(4:6)),bs));
u=bs(ib).v_depart-xp(4:6);
if rand(stream)<.2
 diagnostic.mode='perturbed'; u=u+.08*randn(stream,1,3);
end
xp(4:6)=xp(4:6)+u;
a=ctocscreen.checkArc(xp,finish-tin-w,p);
if ~strcmp(a.status,'ok')||a.min_altitude_km<200, return; end
xa=ctocscreen.propagateTwoBody(xp,da,p.mu_km3_s2);
pool=plan.ids(k+1:min(n,k+4));
events=ctocscreen.scanEncounters(xa,ta,finish-ta,p,pool,cfg.scan_step_s);
events=events(events(:,2)>ta+60 & events(:,2)<finish,:);
if isempty(events), return; end
diagnostic.screen_distance_km=events(1,3);
if events(1,3)>cfg.candidate_gate_km, return; end
diagnostic.gated=1; B=events(1,1); db=events(1,2)-ta;
z0=[w u da db];
limits=[0 u-2 60 60; min(span*.5,oldA*.8) u+2 span*.85 span*.85; finish zeros(1,5)];
diagnostic.refined=1; diagnostic.ids=[A B];
f=ctocscreen.refineFragment(xin,tin,[A B],z0,limits,p,cfg);
if ~f.passed, return; end
% Move discovered B directly behind A and rebuild all affected downstream legs.
j=find(plan.ids==B,1); idx=[1:k j k+1:j-1 j+1:n];
proposal.ids=plan.ids(idx); proposal.reference_v=plan.reference_v(idx,:);
proposal.fixed_dv=plan.fixed_dv(idx,:); proposal.skip=plan.skip(idx);
proposal.waits=plan.waits(idx);
proposal.times(k:k+1)=f.event_times_s;
proposal.waits(k)=f.z(1); proposal.waits(k+1)=0;
proposal.skip(k)=false; proposal.skip(k+1)=true;
proposal.fixed_dv(k,:)=f.z(2:4); proposal.fixed_dv(k+1,:)=NaN;
% Existing S2 pieces after a reordering cannot keep stale fixed impulses.
if j>k+1
 proposal.skip(k+2:end)=false; proposal.fixed_dv(k+2:end,:)=NaN;
end
end
