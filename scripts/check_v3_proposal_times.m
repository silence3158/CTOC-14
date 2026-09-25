function check_v3_proposal_times()
%CHECK_V3_PROPOSAL_TIMES Synthetic proposal boundary and merged-burn checks.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'tests'));
[eph,c,s]=v3SyntheticFixture(); c.joint_target_ids=1;
s.maneuver_times_s=[0;300]; s.delta_v_km_s=[.0001 0 0;0 .0001 0];
p=ctocscreen.v3ShootingProblem(s,eph,c);
trial=p.z0; trial(p.time_indices(2))=-1e-12;
[candidate,report]=p.connect(trial,10);
assert(report.passed&&report.time_projection_max_s>0);
assert(candidate.maneuver_times_s(1)==0);
trial=p.z0; trial(p.time_indices(end))=eph.model.horizon_s/86400+1e-12;
[candidate,report]=p.connect(trial,10);
assert(report.passed&&report.time_projection_max_s>0);
assert(candidate.duration_s==eph.model.horizon_s);
% Coincident proposed impulses are summed and keep the last outgoing anchor.
trial=p.z0; trial(p.time_indices(2))=300/86400;
[candidate,report]=p.connect(trial,10);
assert(report.passed&&report.merged_burns==1);
assert(numel(candidate.maneuver_times_s)==1&&candidate.maneuver_times_s(1)==300);
% A moved encounter uses the nonlinear target position, not its old anchor.
trial=p.z0; epochs=[0;trial(p.time_indices(2:end))*86400];
moved=find(abs(epochs-300)<1e-9); trial(p.time_indices(moved))=315/86400;
[candidate,report]=p.connect(trial,10,'encounter_aligned');
assert(report.passed&&report.aligned_visits==1);
x=ctocscreen.initialState(candidate.initial_q,eph.model.mu,eph.model.re).';
x(4:6)=x(4:6)+candidate.delta_v_km_s(1,:).';
x=ctocscreen.v3Arc(x,0,315,eph.model,c,false,true);
target=ctocscreen.v3QueryTargets(eph,1,315);
assert(norm(x(1:3).'-target)<2*c.connection_tolerance_km);
assert(candidate.witness_times_s(1)==315);
fprintf('PROPOSAL_TIMES_CHECK: boundary retraction and coincident burns passed (synthetic).\n');
end
