function result=test_zhang_dp_tail27(budgetSeconds,label)
%TEST_ZHANG_DP_TAIL27 Authorized historical-order comparison, NOT cold search.
% MATLAB adaptation of Zhang's epoch-pair DP and adaptive tube refinement.
if nargin<1, budgetSeconds=480; end
if nargin<2, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
out=fullfile(sim,'runs/v3/diagnostics',['zhang_dp_tail27_' label]);
assert(~isfolder(out),'Do not overwrite experiment evidence.'); mkdir(out);
diary(fullfile(out,'console.txt')); cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
rng(888,'twister'); timer=tic;
input=fullfile(sim,'runs/v3/diagnostics/tail27_20260925_method_search/comparison.mat');
loaded=load(input,'result'); old=loaded.result;
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim, ...
 'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
m=eph.model; c=ctocscreen.v3Defaults(); c.plane_penalty_km_s=0;
c.budget_s=budgetSeconds; c.guided_max_revolutions=Inf; c.guided_branches=Inf;
sig=ctocscreen.v3ImplementationSignature(eph);
scriptFiles={mfilename('fullpath'),fullfile(sim,'scripts/+zhangdp/solve.m')};
scriptFiles{1}=[scriptFiles{1} '.m'];
scriptHashes=cellfun(@ctocscreen.v3FileHash,scriptFiles,'UniformOutput',false);
order=old.plan(:,1); referenceTimes=old.plan(:,5); original=old.reconstructed_schedule;
firstBurn=original.maneuver_times_s(1);
x0=ctocscreen.initialState(original.initial_q,m.mu,m.re).';
x=ctocscreen.v3Arc(x0,0,firstBurn,m,c,false,true);
result=struct('purpose','user_authorized_teammate_order_diagnostic_not_cold_start', ...
 'upstream_commit','22126db5374f6eb465965a804d1a69706aefe1c0', ...
 'input_file',input,'input_sha256',ctocscreen.v3FileHash(input), ...
 'seed',888,'budget_s',budgetSeconds,'config',c,'source_signature',sig, ...
 'script_files',{scriptFiles},'script_sha256',{scriptHashes},'order',order, ...
 'baseline_verification',old.reconstructed_verification,'baseline_schedule',original, ...
 'trials',{{}},'best_new',[],'elapsed_s',0);
fprintf('AUTHORIZED ORDER COMPARISON baseline %.12f km/s; 35 targets; budget %.1f s\n', ...
 old.reconstructed_verification.total_dv_km_s,budgetSeconds);
% First evaluate the same arrival epochs under the reference code's topology.
centers=referenceTimes; halfWidths=[0 21600 10800 5400 2700 1350 675 337.5];
bestNominal=Inf; reserve=70;
for iteration=1:numel(halfWidths)
 if toc(timer)>=budgetSeconds-reserve, break; end
 width=halfWidths(iteration); trial=struct('iteration',iteration,'half_width_s',width, ...
  'status','started','path',[],'stats',[],'schedule',[],'verification',[], ...
  'elapsed_s',0,'error',''); roundTimer=tic;
 try
  times=cell(36,1); positions=cell(36,1); times{1}=firstBurn; positions{1}=x(1:3).';
  for k=1:35
   if width==0, values=centers(k); else, values=centers(k)+linspace(-width,width,9); end
   values=unique(values(values>firstBurn+1 & values<=m.horizon_s));
   times{k+1}=values; positions{k+1}=ctocscreen.v3QueryTargets(eph,order(k),values,'pairs');
  end
  [path,stats]=zhangdp.solve(times,positions,x(4:6).',m.mu,m.re,@()toc(timer)>budgetSeconds-reserve);
  trial.path=path; trial.stats=stats; trial.status='proposal_only';
  fprintf('DP iteration %d width %.3f s proposal %.9f km/s pairs %d branches %d time %.3f s\n', ...
   iteration,width,path.nominal_dv_km_s,stats.pairs,stats.branches,stats.elapsed_s);
  if path.nominal_dv_km_s<bestNominal
   bestNominal=path.nominal_dv_km_s; centers=path.arrival_times_s;
   % Preserve each improving graph result before potentially failing J2 work.
   result.pending_path=path; save(fullfile(out,'experiment.mat'),'result','-v7.3');
   [schedule,detail]=restore(path,original.initial_q,order,x,eph,c,timer);
   trial.schedule=schedule; trial.repair_detail=detail;
   [verification,~]=ctocscreen.v3Replay(schedule,eph,c,true);
   trial.verification=verification; trial.status=verification.status;
   fprintf('INDEPENDENT iteration %d passed %d visits %d dv %.12f height %.6f distance %.9f elapsed %.3f\n', ...
    iteration,verification.passed,verification.visit_count,verification.total_dv_km_s, ...
    verification.min_altitude_lower_km,max(verification.distance_km),toc(timer));
   if verification.passed && (isempty(result.best_new)|| ...
     verification.total_dv_km_s<result.best_new.verification.total_dv_km_s)
    result.best_new=struct('schedule',schedule,'verification',verification,'iteration',iteration);
    save(fullfile(out,'best_new.mat'),'-struct','result','best_new');
   end
  end
 catch err
  trial.status='failed'; trial.error=[err.identifier ': ' err.message];
  fprintf('FAIL iteration %d %s\n',iteration,trial.error);
 end
 trial.elapsed_s=toc(roundTimer); result.trials{end+1}=trial;
 result.elapsed_s=toc(timer); save(fullfile(out,'experiment.mat'),'result','-v7.3');
end
assert(isequal(sig,ctocscreen.v3ImplementationSignature(eph)),'Sources changed during experiment.');
assert(isequal(scriptHashes,cellfun(@ctocscreen.v3FileHash,scriptFiles,'UniformOutput',false)), ...
 'Experiment sources changed while running.');
result.elapsed_s=toc(timer); result.source_unchanged=true;
if ~isempty(result.best_new)
 result.improvement_km_s=old.reconstructed_verification.total_dv_km_s-result.best_new.verification.total_dv_km_s;
 fprintf('FINAL new %.12f km/s improvement %.12f km/s elapsed %.3f s\n', ...
  result.best_new.verification.total_dv_km_s,result.improvement_km_s,result.elapsed_s);
else
 fprintf('FINAL no independently verified new trajectory elapsed %.3f s\n',result.elapsed_s);
end
save(fullfile(out,'experiment.mat'),'result','-v7.3');
fprintf('SAVED %s\n',out);
end

function [s,detail]=restore(path,q,order,x,eph,c,timer)
s=struct('schema_version','free_maneuver_v3','dynamics_id','central_j2', ...
 'initial_q',q,'duration_s',path.arrival_times_s(end), ...
 'maneuver_times_s',path.departure_times_s,'delta_v_km_s',zeros(35,3), ...
 'witness_times_s',nan(35,1));
s.witness_times_s(order)=path.arrival_times_s; detail=cell(35,1);
for k=1:35
 if toc(timer)>c.budget_s-15, error('zhangdp:repairBudget','Repair budget exhausted.'); end
 goal=ctocscreen.v3QueryTargets(eph,order(k),path.arrival_times_s(k),'pairs');
 [dv,detail{k}]=ctocscreen.v3GuidedTransfer(x,path.departure_times_s(k), ...
  path.arrival_times_s(k),goal,eph.model,c,timer,path.v_depart(k,:));
 s.delta_v_km_s(k,:)=dv.';
 x=ctocscreen.v3Arc([x(1:3);x(4:6)+dv],path.departure_times_s(k), ...
  path.arrival_times_s(k),eph.model,c,false,true);
end
end
