function diagnosis=diagnose_v4_first(label)
%DIAGNOSE_V4_FIRST Read failed controls only to identify numerical defects.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
path=fullfile(sim,'runs/v4/development',['first_failure_diagnosis_' label '.mat']);
assert(~isfile(path),'Diagnosis output already exists.');
s=load(fullfile(sim,'runs/v4/search/v4_cold_first_300_seed888_01/result.mat'));
r=s.result; c=r.manifest.config;
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
q=r.best.q;
[a,tr]=ctocscreen.v4.replay(q,eph,c);
[v,vr]=ctocscreen.v4.replay(q,eph,c,true);
assert(~strcmp(a.status,'propagation_failure'),a.failure_reason);
assert(~strcmp(v.status,'propagation_failure'),v.failure_reason);
fprintf('Saved screen %d; fresh screen %d; independent %d; cached/fresh endpoint %.9g km\n', ...
 r.best.actual.visit_count,a.visit_count,v.visit_count,norm(a.final_state(1:3)-r.best.actual.final_state(1:3)));
rows=zeros(numel(tr.arcs),4);
for j=1:numel(tr.arcs)
 ta=tr.arc_start(j); tb=tr.arc_end(j);
 x=ctocscreen.v4.stateAt(tr,ta,'post');
 y=ctocscreen.v3Arc(x,ta,tb,eph.model,c,false,true);
 yscreen=ctocscreen.v4.stateAt(tr,tb,'pre'); yref=ctocscreen.v4.stateAt(vr,tb,'pre');
 rows(j,:)=[tb,norm(y(1:3)-yscreen(1:3)),norm(yref(1:3)-yscreen(1:3)),norm(yref(4:6)-yscreen(4:6))];
end
disp(array2table(rows,'VariableNames',{'time_s','one_arc_error_km','chain_error_km','chain_error_km_s'}));
node=r.best;
% An early physical prefix distinguishes scaling/event bugs from long-chain sensitivity.
node.q.T=q.tau(6); node.q.tau=q.tau(1:5); node.q.u=q.u(1:5,:);
node.q.witness(q.witness>node.q.T)=NaN;
[node.actual,node.trace]=ctocscreen.v4.replay(node.q,eph,c);
ids=find(node.actual.distance_km<=c.search_radius_km);
p=ctocscreen.v4.problem(node,ids,node.actual.witness_times_s(ids),'full',node.q.T,eph,c);
z=p.z0; e=ctocscreen.v4.evaluate(p,z,true); steps=cell(1,7);
fprintf('Prefix visits %d; J %.9g; V %.9g; n=%d; maxJac %.9g\n',numel(ids),e.J,e.V,p.n,max(abs(nonzeros(e.Jeq))));
for j=1:7
 trust=.03/4^(j-1); st=ctocscreen.v4.coneStep(p,z,e,trust,10,true,3);
 st.trust=trust;
 if ~isempty(st.d)
  [pred,st.predicted_V,st.predicted_J]=ctocscreen.v4.modelMerit(p,z,e,st.d,10,true);
  st.predicted_drop=e.J/p.sv+10*e.V-pred;
  st.order_residual=max(p.A*(z+st.d));
  try
   trial=ctocscreen.v4.evaluate(p,z+st.d,false);
   st.actual_V=trial.V; st.actual_drop=e.J/p.sv+10*e.V-(trial.J/p.sv+10*trial.V);
   st.ratio=st.actual_drop/st.predicted_drop;
  catch err
   st.failure_id=err.identifier; st.failure_reason=err.message;
  end
 end
 disp(rmfield(st,{'d','output'})); steps{j}=st;
end
diagnosis=struct('diagnostic_only',true,'source_signature',ctocscreen.v4.signature(), ...
 'screen',a,'independent',v,'arc_comparison',rows,'steps',{steps});
save(path,'diagnosis');
end
