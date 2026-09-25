function report=report_v3_cost_progress(label)
%REPORT_V3_COST_PROGRESS Diagnose saved runs only; never supply search inputs.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
folder=fullfile(sim,'runs','v3','search',label);
d=load(fullfile(folder,'checkpoint.mat'),'state'); state=d.state;
report=report_v3_cold_run(label);
fprintf('Recovery=%d replanQueue=%d replans=%d pendingFresh=%d\n', ...
 numel(state.recovery_pool),numel(state.replan_queue),numel(state.replan_history),numel(state.fresh_queue));
for k=1:numel(state.replan_history)
 r=state.replan_history{k};
 fprintf('Replan %d %s cuts=%s elapsed=%.3f\n',k,r.status,mat2str(r.cut_times_s),r.elapsed_s);
 for j=1:numel(r.attempts)
  a=r.attempts{j}; spent=NaN;
  if isfield(a,'retained_dv_km_s'), spent=a.retained_dv_km_s; end
  fprintf('  %s retained=%.6f proposals=%d %s\n',a.mode,spent,a.proposals,a.reason);
 end
end
if ~isempty(state.best_rejected_complete), s=state.best_rejected_complete.schedule;
elseif ~isempty(state.best_partial), s=state.best_partial.schedule;
else, return; end
dv=vecnorm(s.delta_v_km_s,2,2); cumulative=cumsum(dv);
fprintf('Saved diagnostic trajectory, NOT initialization:\n');
for k=1:numel(dv)
 ta=s.maneuver_times_s(k); tb=s.duration_s;
 if k<numel(dv), tb=s.maneuver_times_s(k+1); end
 ids=find(s.witness_times_s>ta&s.witness_times_s<=tb);
 fprintf('  burn%02d t=%8.1f dt=%8.1f dv=%8.5f sum=%8.5f targets=%s\n', ...
  k,ta,tb-ta,dv(k),cumulative(k),mat2str(ids.'));
end
for k=1:numel(state.history)
 h=state.history{k};
 fprintf('Joint%02d %s N=%d status=%s selected=%s error=%s\n',k,h.operator, ...
  numel(h.task_target_ids),h.diagnostic.status,h.diagnostic.selected_source,h.diagnostic.failure_reason);
 info=h.diagnostic;
 fuelStarted=isfield(info,'fuel_started')&&info.fuel_started;
 fprintf('  elapsed=%.3f work=%d fuelStarted=%d\n',info.elapsed_s,info.refinement_performed,fuelStarted);
 if isfield(info,'accepted_steps'), fprintf('  acceptedSteps=%d trials=%d\n',info.accepted_steps,numel(info.trials)); end
 if isfield(info,'trials')
  for j=1:numel(info.trials)
   trial=info.trials{j};
   fprintf('  trial%d connected=%d accepted=%d proposalDV=%.8g %s\n', ...
    j,trial.connection.passed,trial.accepted,trial.actual_dv_km_s,trial.connection.reason);
   if isfield(trial.connection,'time_projection_max_s')
    fprintf('    timeProjection=%.9g minRawGap=%.9g mode=%s aligned=%d\n', ...
     trial.connection.time_projection_max_s,trial.connection.raw_min_event_gap_s, ...
     trial.connection.projection_mode,trial.connection.aligned_visits);
   end
   if ~isempty(trial.verification)
    fprintf('    verification=%s activeMiss=%.8g\n',trial.verification.status, ...
     max(trial.verification.distance_km(info.active_target_ids)));
   end
  end
 end
 if isfield(info,'fuel')&&isfield(info.fuel,'exitflag'), fprintf('  fuelExitflag=%g\n',info.fuel.exitflag); end
end
if ~isempty(state.failures)
 ids=cellfun(@(f)f.id,state.failures,'UniformOutput',false); [names,~,index]=unique(ids);
 for k=1:numel(names), fprintf('Failure %d %s\n',sum(index==k),names{k}); end
end
end
