function result=test_zhang_delayed_tail27(budgetSeconds,label)
%TEST_ZHANG_DELAYED_TAIL27 Explicit-wait adaptation; historical diagnostic only.
if nargin<1, budgetSeconds=300; end
if nargin<2, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
out=fullfile(sim,'runs/v3/diagnostics',['zhang_delayed_tail27_' label]);
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
files={[mfilename('fullpath') '.m'],fullfile(sim,'scripts/+zhangdp/delayed.m')};
hashes=cellfun(@ctocscreen.v3FileHash,files,'UniformOutput',false);
original=old.reconstructed_schedule; order=old.plan(:,1); centers=old.plan(:,5); waits=old.plan(:,2);
first=original.maneuver_times_s(1); x0=ctocscreen.initialState(original.initial_q,m.mu,m.re).';
x=ctocscreen.v3Arc(x0,0,first,m,c,false,true);
result=struct('purpose','authorized_historical_order_delayed_comparison_not_cold_start', ...
 'upstream_commit','22126db5374f6eb465965a804d1a69706aefe1c0', ...
 'input_file',input,'input_sha256',ctocscreen.v3FileHash(input),'seed',888, ...
 'budget_s',budgetSeconds,'config',c,'source_signature',sig,'script_files',{files}, ...
 'script_sha256',{hashes},'baseline_verification',old.reconstructed_verification, ...
 'baseline_schedule',original,'order',order,'trials',{{}},'best_new',[],'elapsed_s',0);
fprintf('DELAYED baseline %.12f km/s budget %.1f s\n',old.reconstructed_verification.total_dv_km_s,budgetSeconds);
widths=[0 21600 10800 5400 2700 1350 675 337.5];
waitWidths=[0 5400 2700 1350 675 337.5 168.75 84.375]; bestNominal=Inf; incumbent=[];
for it=1:numel(widths)
 if toc(timer)>=budgetSeconds-65, break; end
 trial=struct('iteration',it,'width_s',widths(it),'wait_width_s',waitWidths(it), ...
  'status','started','path',[],'stats',[],'schedule',[],'verification',[],'error','');
 try
  times=cell(35,1);
  for k=1:35
   if widths(it)==0, t=centers(k); else, t=centers(k)+linspace(-widths(it),widths(it),5); end
   times{k}=unique(t(t>first+1&t<=m.horizon_s));
  end
  [path,stats]=zhangdp.delayed(times,order,eph,x,first,waits,waitWidths(it),12, ...
   @()toc(timer)>budgetSeconds-65,incumbent);
  trial.path=path; trial.stats=stats; trial.status='proposal_only';
  fprintf('DELAYED DP iteration %d proposal %.12f pairs %d branches %d elapsed %.3f\n', ...
   it,path.nominal_dv_km_s,stats.pairs,stats.branches,toc(timer));
  if path.nominal_dv_km_s<bestNominal
   bestNominal=path.nominal_dv_km_s; centers=path.arrival_times_s; incumbent=path;
   waits=[first;path.departure_times_s(2:end)-path.arrival_times_s(1:end-1)];
   result.pending_path=path; save(fullfile(out,'experiment.mat'),'result','-v7.3');
   [s,detail]=repair(path,original.initial_q,order,x,eph,c,timer);
   trial.schedule=s; trial.repair_detail=detail;
   [v,~]=ctocscreen.v3Replay(s,eph,c,true); trial.verification=v; trial.status=v.status;
   fprintf('DELAYED INDEPENDENT iteration %d passed %d visits %d dv %.12f height %.6f distance %.9f elapsed %.3f\n', ...
    it,v.passed,v.visit_count,v.total_dv_km_s,v.min_altitude_lower_km,max(v.distance_km),toc(timer));
   if v.passed && (isempty(result.best_new)||v.total_dv_km_s<result.best_new.verification.total_dv_km_s)
    result.best_new=struct('schedule',s,'verification',v,'iteration',it);
    save(fullfile(out,'best_new.mat'),'-struct','result','best_new');
   end
  end
 catch err
  trial.status='failed'; trial.error=[err.identifier ': ' err.message];
  fprintf('DELAYED FAIL iteration %d %s\n',it,trial.error);
 end
 result.trials{end+1}=trial; result.elapsed_s=toc(timer); save(fullfile(out,'experiment.mat'),'result','-v7.3');
end
assert(isequal(sig,ctocscreen.v3ImplementationSignature(eph)),'Sources changed during run.');
assert(isequal(hashes,cellfun(@ctocscreen.v3FileHash,files,'UniformOutput',false)), ...
 'Experiment sources changed during run.');
result.elapsed_s=toc(timer); result.source_unchanged=true;
if ~isempty(result.best_new)
 result.improvement_km_s=old.reconstructed_verification.total_dv_km_s-result.best_new.verification.total_dv_km_s;
 fprintf('DELAYED FINAL new %.12f improvement %.12f elapsed %.3f\n', ...
  result.best_new.verification.total_dv_km_s,result.improvement_km_s,result.elapsed_s);
else
 fprintf('DELAYED FINAL no independently verified new trajectory elapsed %.3f\n',result.elapsed_s);
end
save(fullfile(out,'experiment.mat'),'result','-v7.3'); fprintf('SAVED %s\n',out);
end

function [s,detail]=repair(path,q,order,x,eph,c,timer)
s=struct('schema_version','free_maneuver_v3','dynamics_id','central_j2', ...
 'initial_q',q,'duration_s',path.arrival_times_s(end), ...
 'maneuver_times_s',path.departure_times_s,'delta_v_km_s',zeros(35,3), ...
 'witness_times_s',nan(35,1));
s.witness_times_s(order)=path.arrival_times_s; detail=cell(35,1); t=path.departure_times_s(1);
for k=1:35
 if toc(timer)>c.budget_s-15, error('zhangdp:repairBudget','Repair budget exhausted.'); end
 td=path.departure_times_s(k); ta=path.arrival_times_s(k);
 if td>t, x=ctocscreen.v3Arc(x,t,td,eph.model,c,false,true); end
 goal=ctocscreen.v3QueryTargets(eph,order(k),ta,'pairs');
 [dv,detail{k}]=ctocscreen.v3GuidedTransfer(x,td,ta,goal,eph.model,c,timer,path.v_depart(k,:));
 s.delta_v_km_s(k,:)=dv.';
 x=ctocscreen.v3Arc([x(1:3);x(4:6)+dv],td,ta,eph.model,c,false,true); t=ta;
end
end
