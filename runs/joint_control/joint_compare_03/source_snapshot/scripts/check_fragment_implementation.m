function check_fragment_implementation()
root=fileparts(fileparts(mfilename('fullpath'))); cd(root); addpath('src');
res=runtests('tests/TestFragments.m'); disp(table(res)); assert(all([res.Passed]));
src=load('runs/screening/campaign16_20260921/elite.mat');
s=ctocscreen.importV1Candidate(src.elite.candidate,src.elite.evaluation);
r=ctocscreen.propagateSchedule(s,src.problem,true); assert(r.passed,r.failure_reason);
assert(abs(r.total_dv_km_s-src.elite.total_dv_km_s)<1e-10);
fprintf('V1_IMPORT_PASS dv%.12f visits%d error_m%.8f min_alt%.6f\n', ...
 r.total_dv_km_s,r.unique_visit_count,1000*max(r.event_distances_km),r.min_altitude_km);
if ~isfolder('runs/fragments/validation'), mkdir('runs/fragments/validation'); end
save('runs/fragments/validation/checks.mat','res','r','s');
end
