function analyze_v3_budget_runs()
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
for label={'v3_budget_480_01','v3_budget_960_01'}
 a=load(fullfile(sim,'runs/v3/search',label{1},'checkpoint.mat')); s=a.state; c=s.config;
 fprintf('\nRUN %s elapsed %.3f iter %d seed %d joint %.1f prefix %.1f restoration %d fueliter %d\n',label{1},s.elapsed_s,s.iteration,c.seed,c.joint_seconds,c.prefix_joint_seconds,c.restoration_iterations,c.joint_iterations);
 fprintf('signature_current %d roots %d work %d elite %d failures %d\n',isequaln(s.signature,ctocscreen.v3ImplementationSignature(eph)),numel(s.root_archive),numel(s.work_pool),~isempty(s.elite),numel(s.failures));
 constructionTime=0; trunc=0; prune=0;
 for k=1:numel(s.construction_history)
  h=s.construction_history{k}; constructionTime=constructionTime+h.elapsed_s;
  trunc=trunc+h.consistency.truncated; prune=prune+h.pruned_extensions;
  fprintf('C%02d t %.1f visits %d>%d raw %d trunc %d prune %d\n',k,h.started_at_s,h.start_visits,h.end_visits,h.incremental_visits,h.consistency.truncated,h.pruned_extensions);
 end
 jointTime=0;
 for k=1:numel(s.history)
  h=s.history{k}; d=h.diagnostic; jointTime=jointTime+d.elapsed_s;
  fprintf('H%02d iter %d %s %s source %s sec %.2f fuel %d flag %g',k,h.iteration,h.operator,d.status,d.selected_source,d.elapsed_s,d.fuel_started,d.exitflag);
  if isfield(d,'initial_violation'), fprintf(' violation %.4g>%.4g',d.initial_violation,d.max_scaled_violation); end
  if ~isempty(d.replay), fprintf(' count %d J %.5f',d.replay.visit_count,d.replay.total_dv_km_s); end
  fprintf('\n');
  if ~isempty(d.restoration), fprintf(' restore iterations %d merit %.4g>%.4g\n',d.restoration.output.iterations,d.restoration.initial_merit,d.restoration.final_merit); end
  if isfield(d,'solver_output'), disp(d.solver_output); end
  if ~isempty(d.failure_reason), fprintf('ERROR %s\n',d.failure_reason); end
 end
 fprintf('TIMES construct %.3f joint %.3f other %.3f trunc %d/%d prune %d\n',constructionTime,jointTime,s.elapsed_s-constructionTime-jointTime,trunc,numel(s.construction_history),prune);
 b=s.best_partial; r=ctocscreen.v3Replay(b.schedule,eph,c,true);
 fprintf('BEST count %d independent %d height %d J %.9f duration %.3f burns %d\n',sum(b.visited),r.visit_count,r.height_passed,b.J,b.t,numel(b.schedule.maneuver_times_s));
 disp(b.schedule.initial_q);
 ids=cellfun(@(f)f.id,s.failures,'UniformOutput',false); [u,~,ix]=unique(ids);
 for j=1:numel(u), fprintf('FAIL %s %d\n',u{j},sum(ix==j)); end
 save(fullfile(sim,'runs/v3/development',[label{1} '_analysis.mat']),'r','constructionTime','jointTime','trunc','prune');
end
end
