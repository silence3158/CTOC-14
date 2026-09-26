function result=report_v4_run(folder)
%REPORT_V4_RUN Read actual independent results, not plan length or penalized cost.
if nargin<1
 sim=fileparts(fileparts(mfilename('fullpath'))); folders=dir(fullfile(sim,'runs','v4','search','*'));
 folders=folders([folders.isdir]&~ismember({folders.name},{'.','..'}));
 assert(~isempty(folders),'No V4 run exists.'); [~,j]=max([folders.datenum]); folder=fullfile(folders(j).folder,folders(j).name);
end
data=load(fullfile(folder,'result.mat'),'result'); result=data.result; a=result.verification;
fprintf('Folder: %s\nCold start: %d; seed: %d\n',folder,result.manifest.cold_start,result.manifest.config.seed);
fprintf('Independent: %d/35; raw delta-V: %.12f km/s; complete pass: %d\n',a.visit_count,a.total_dv_km_s,a.passed);
fprintf('Height lower bound: %.9f km; initial orbit: %d; time: %.3f s\n',a.min_altitude_lower_km,a.initial_passed,result.best.q.T);
fprintf('Actual elapsed: %.3f s; source unchanged: %d; ATK aligned: %d\n',result.stats.total_seconds,result.stats.source_unchanged,a.official_alignment_verified);
fprintf('B calls shared/full/structure: %d/%d/%d; total iterations: %d\n',result.stats.shared_calls,result.stats.full_calls,result.stats.structure_calls,result.stats.joint_iterations);
end
