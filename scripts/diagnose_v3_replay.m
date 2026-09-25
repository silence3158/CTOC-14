function diagnose_v3_replay()
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
a=load(fullfile(sim,'runs/v3/search/v3_user_01_revised/checkpoint.mat')); c=ctocscreen.v3Defaults(a.state.config);
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
s=a.state.best_partial.schedule; m=eph.model;
[r,tr]=ctocscreen.v3Replay(s,eph,c,false);
[v,tv]=ctocscreen.v3Replay(s,eph,c,true);
fprintf('screen %d independent %d\n',r.visit_count,v.visit_count);
for k=1:numel(tr.arcs)
 so=tr.arcs{k}; si=tv.arcs{k}; xs=deval(so,so.x(end)); xi=deval(si,si.x(end));
 fprintf('arc %d t %.3f dt %.3f dr %.9g dv %.9g R %.6g\n',k,so.x(end),diff(so.x([1 end])),norm(xs(1:3)-xi(1:3)),norm(xs(4:6)-xi(4:6)),norm(xs(1:3)));
end
p=ctocscreen.v3ShootingProblem(s,eph,c); [ci,eq]=p.constraints(p.z0);
fprintf('NLP initial violation %.9g expanded bounds %d\n',max([0;ci;abs(eq)]),p.expanded_bound_count);
for j=[1 10 20 35]
 t=s.witness_times_s(j); k=find(cellfun(@(a)t>=a.x(1)&&t<=a.x(end),tr.arcs),1); x=deval(tr.arcs{k},t);
 xn=p.z0(p.state_indices(:,p.witness_nodes(j))).*[1e4;1e4;1e4;1;1;1];
 fprintf('witness %d t %.9g node dr %.9g\n',j,t,norm(x(1:3)-xn(1:3)));
end
[stable,sr,stability]=ctocscreen.v3StablePrefix(s,eph,c);
fprintf('Stable prefix visits %d duration %.9g J %.9g\n',sr.visit_count,stable.duration_s,sr.total_dv_km_s);
disp(stability);
save(fullfile(sim,'runs/v3/development/replay_diagnosis.mat'),'r','v','tr','tv','ci','eq','stable','sr','stability');
end
