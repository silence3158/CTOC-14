function check_v3_connected_mutation(label)
% Bounded development check of a different physical initial orbit.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
if nargin<1, label='core_connected_mutation'; end
assert(~isempty(regexp(label,'^[A-Za-z0-9_-]+$','once')));
w=load(fullfile(sim,'runs/v3/search/warm_repair_01/warm_start.mat'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v3Defaults(struct('shooting_reltol',1e-13,'max_step_s',60));
[s,name,connection]=ctocscreen.v3ConnectedMutation(w.s,eph,c,RandStream('mt19937ar','Seed',101),30,4);
verification=[];
if connection.passed, verification=ctocscreen.v3Verify(s,eph,c); end
save(fullfile(sim,'runs/v3/development',[label '.mat']),'s','name','connection','verification','c');
disp(connection); fprintf('initial q change: '); fprintf('%.6g ',s.initial_q-w.s.initial_q); fprintf('\n');
if ~isempty(verification)
 fprintf('MUTATION J %.12f visits %d passed %d maxdist %.9g\n',verification.total_dv_km_s, ...
  verification.visit_count,verification.passed,max(verification.distance_km));
end
end
