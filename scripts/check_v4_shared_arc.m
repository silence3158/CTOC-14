function result=check_v4_shared_arc(label)
%CHECK_V4_SHARED_ARC One decision-bearing component check before warm search.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
if nargin<1, label='02'; end
folder=fullfile(sim,'runs/v4/development',['shared_arc_integration_check_' label]);
assert(~isfolder(folder)); mkdir(folder); clock=tic;
c=ctocscreen.v4.defaults(); sig=ctocscreen.v4.signature();
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
data=load(fullfile(sim,'runs/v4/development/multiarc_probe_minimal_03/report.mat'));
s=data.report.records{7}.sample;
[a,tr]=ctocscreen.v4.replay(s.q,eph,c,false);
node=struct('q',tr.q,'actual',a,'trace',tr,'generation',numel(tr.q.tau), ...
 'attempts',0,'zero_gain',0,'heuristic_H',NaN,'origin','diagnostic'); node.q.witness=a.witness_times_s;
[a,tr]=ctocscreen.v4.replay(s.parent,eph,c,false); parent=node;
parent.q=tr.q; parent.q.witness=a.witness_times_s; parent.actual=a; parent.trace=tr; parent.generation=node.generation-1;
[children,report]=ctocscreen.v4.sharedArc(node,parent,eph,c,7);
assert(~isempty(children),'Known missed encounter was not recovered.');
result=struct('report',report,'verified',{{}},'signature',sig,'input',s);
for k=1:numel(children)
 child=children{k}; [v,~]=ctocscreen.v4.replay(child.q,eph,c,true);
 active=node.actual.distance_km<=1;
 assert(v.independent&&v.initial_passed&&v.height_passed&&all(v.distance_km(active)<=1));
 fprintf('CHILD %d visits=%d added=%s target22=%.6f J=%.9f\n',k,v.visit_count,mat2str(find(v.distance_km<=1&~active).'),v.distance_km(22),v.total_dv_km_s);
 assert(v.visit_count>=21&&numel(child.q.tau)==numel(node.q.tau));
 M=numel(parent.q.tau);
 assert(isequal(child.q.u(1:M,:),parent.q.u)&&isequal(child.q.tau(1:M),parent.q.tau));
 assert(child.q.tau(end)>=parent.q.T&&child.generation==parent.generation+1);
 result.verified{end+1}=struct('q',child.q,'verification',v);
end
result.seconds=toc(clock); result.source_unchanged=isequal(sig,ctocscreen.v4.signature()); assert(result.source_unchanged);
save(fullfile(folder,'result.mat'),'result');
assert(any(cellfun(@(s)s.verification.distance_km(22)<=1,result.verified)),'Target22 not recovered among alternatives.');
fprintf('SHARED CHECK children=%d attempts=%d far_accepted=%d seconds=%.3f\n',numel(children),report.attempts,report.far_accepted,result.seconds);
end
