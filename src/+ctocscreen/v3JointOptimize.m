function [candidate,info]=v3JointOptimize(seed,eph,c)
%V3JOINTOPTIMIZE Joint feasibility restoration, then fuel, same variables.
started=tic; candidate=seed;
info=struct('status','failed','exitflag',NaN,'failure_reason','', ...
 'max_scaled_violation',Inf,'elapsed_s',0,'replay',[],'selected_source','solver_candidate', ...
 'active_target_ids',c.joint_target_ids,'restoration',[],'fuel_started',false,'refinement_performed',false);
baseline=[]; continuation=[];
try
 independentSeed=all(isfinite(seed.witness_times_s));
 [baseline,baselineTrace]=ctocscreen.v3JointReplay(seed,eph,c,independentSeed);
 assert(~strcmp(baseline.status,'propagation_failure'),'ctocscreen:v3:seedReplay','Seed propagation failed.');
 info.seed_visit_count=baseline.visit_count;
 if c.projected_joint && ctocscreen.v3RecoveryMerit(baseline,c.joint_target_ids,1)==0 && ...
    ~isempty(seed.maneuver_times_s) && ~isfield(seed,'visit_plan_times_s')
  cached=[];
  if independentSeed||c.shooting_reltol<=1e-13
   cached=struct('schedule_key',ctocscreen.v3ScheduleKey(seed),'replay',baseline,'trace',baselineTrace);
  end
  [candidate,info]=ctocscreen.v3ProjectedJoint(seed,eph,c,baseline,max(0,c.joint_seconds-toc(started)),cached);
  info.elapsed_s=toc(started); info.budget_overrun_s=max(0,info.elapsed_s-c.joint_seconds);
  return
 end
 joint=c; joint.shooting_reltol=min(c.shooting_reltol,1e-12);
 p=ctocscreen.v3ShootingProblem(seed,eph,joint);
 context=struct('ephemeris_signature',eph.signature,'model',eph.model,'states0',eph.states0, ...
  'search_radius_km',c.search_radius_km,'shooting_reltol',joint.shooting_reltol,'max_step_s',c.max_step_s);
 start=p.z0; resumed=false;
 if isfield(seed,'restoration_state')
  saved=seed.restoration_state;
  if isfield(saved,'context')&&isequaln(saved.context,context) ...
    &&isequal(saved.target_ids,c.joint_target_ids)&&isequal(saved.source_z,p.z0) ...
    &&numel(saved.iterate)==p.variable_count&&all(isfinite(saved.iterate)) ...
    &&all(saved.iterate>=p.lb)&&all(saved.iterate<=p.ub)
   start=saved.iterate; resumed=true;
  end
 end
 [c0,e0]=p.constraints(p.z0); assert(all(isfinite([c0;e0])));
 info.initial_violation=max([0;c0;abs(e0)]); info.variable_count=p.variable_count;
 info.expanded_bound_count=p.expanded_bound_count;
 info.initial_defects=defects(e0); z=p.z0;
 [z,info.restoration]=ctocscreen.v3Restore(p,c,max(0,c.joint_seconds*c.restoration_fraction-toc(started)),start);
 info.restoration.resumed=resumed;
 continuation=struct('target_ids',c.joint_target_ids,'source_z',p.z0,'iterate',z,'context',context);
 info.refinement_performed=info.restoration.refinement_performed;
 info.exitflag=info.restoration.exitflag; info.solver_output=info.restoration.output;
 [cc,ee]=p.constraints(z); violation=max([0;cc;abs(ee);p.A*z-p.b]);
 info.restoration.final_defects=defects(ee);
 if violation<=10*c.constraint_tolerance && toc(started)<c.joint_seconds
  info.fuel_started=true;
  beforeFuel=z;
  [zf,raw,info.fuel]=ctocscreen.v3Fuel(p,z,c,max(0,.7*(c.joint_seconds-toc(started))));
  info.refinement_performed=info.refinement_performed||info.fuel.refinement_performed;
  flag=info.fuel.exitflag; output=info.fuel.output;
  [cf,ef]=p.constraints(zf); vf=max([0;cf;abs(ef);p.A*zf-p.b]);
  if vf<=10*c.constraint_tolerance, z=zf; cc=cf; ee=ef; end
  info.exitflag=flag; info.solver_output=output;
  % A multiple-shooting step is only a proposal. Backtrack its physical
  % controls under fixed-pulse independent replay, never reset node states.
  info.fuel.physical_trials=0; info.fuel.physical_accepted=false;
  if p.objective(raw)<p.objective(beforeFuel)
   if baseline.passed, baseIndependent=baseline;
   else, baseIndependent=ctocscreen.v3JointReplay(p.decode(beforeFuel),eph,c,true); end
   for alpha=[1 .3 .1 .03 .01 .003 .001]
    if toc(started)>=c.joint_seconds, break; end
    trial=beforeFuel+alpha*(raw-beforeFuel);
    check=ctocscreen.v3JointReplay(p.decode(trial),eph,c,true); info.fuel.physical_trials=info.fuel.physical_trials+1;
    if ctocscreen.v3RecoveryMerit(check,c.joint_target_ids,1)==0 ...
       &&check.visit_count>=baseIndependent.visit_count&&check.total_dv_km_s<baseIndependent.total_dv_km_s
     z=trial; [cc,ee]=p.constraints(z); info.fuel.physical_accepted=true;
     info.fuel.independent=check; info.fuel.alpha=alpha; break
    end
   end
  end
 end
 info.max_scaled_violation=max([0;cc;abs(ee);p.A*z-p.b]); info.final_defects=defects(ee);
 candidate=ctocscreen.v3Normalize(p.decode(z),eph.model);
 info.replay=ctocscreen.v3JointReplay(candidate,eph,c,all(isfinite(seed.witness_times_s)));
 candidate.witness_times_s=info.replay.witness_times_s;
 candidate.witness_times_s(info.replay.distance_km>1)=NaN;
 candidate.validation_level='proposal'; info.status='not_feasible';
 if info.replay.passed
  info.status='screened_feasible'; candidate.validation_level='nominal_j2_screened';
 elseif ctocscreen.v3RecoveryMerit(info.replay,c.joint_target_ids,1)==0
  info.status='prefix_restored';
 end
catch err
 info.failure_reason=err.message; info.failure_id=err.identifier;
end
% Infeasible recovery is ordered by physical miss, not by fuel or visit count.
before=ctocscreen.v3RecoveryMerit(baseline,c.joint_target_ids,c.search_radius_km);
after=ctocscreen.v3RecoveryMerit(info.replay,c.joint_target_ids,c.search_radius_km);
baselineFeasible=ctocscreen.v3RecoveryMerit(baseline,c.joint_target_ids,1)==0;
proposalFeasible=ctocscreen.v3RecoveryMerit(info.replay,c.joint_target_ids,1)==0;
retain=before<=after;
if baselineFeasible
 retain=~proposalFeasible||baseline.total_dv_km_s<=info.replay.total_dv_km_s;
elseif proposalFeasible
 retain=false;
end
if ~isempty(baseline)&&baseline.height_passed&&retain
 candidate=seed; candidate.witness_times_s=baseline.witness_times_s;
 candidate.witness_times_s(baseline.distance_km>1)=NaN;
 info.replay=baseline; info.selected_source='retained_recovery_seed';
 if baseline.passed
  candidate.validation_level='nominal_j2_screened'; info.status='screened_feasible';
  info.selected_source='retained_feasible_seed';
 elseif baselineFeasible
  info.status='prefix_restored'; info.selected_source='retained_prefix';
 end
end
if isfield(candidate,'restoration_state'), candidate=rmfield(candidate,'restoration_state'); end
info.retained_auxiliary_state=false;
if numel(c.joint_target_ids)==35&&~isempty(continuation)&&strcmp(info.selected_source,'retained_recovery_seed')
 candidate.restoration_state=continuation; info.retained_auxiliary_state=true;
end
info.recovery_initial_merit=before;
info.recovery_final_merit=ctocscreen.v3RecoveryMerit(info.replay,c.joint_target_ids,c.search_radius_km);
info.recovery_progress=info.recovery_final_merit<before;
info.elapsed_s=toc(started); info.budget_overrun_s=max(0,info.elapsed_s-c.joint_seconds);
end

function d=defects(eq)
x=reshape(eq,6,[]); d=struct('max_position_km',max(vecnorm(x(1:3,:)*1e4,2,1)), ...
 'max_velocity_km_s',max(vecnorm(x(4:6,:),2,1)));
end
