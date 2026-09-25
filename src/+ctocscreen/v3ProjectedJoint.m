function [candidate,info]=v3ProjectedJoint(seed,eph,c,baseline,budget,cached)
%V3PROJECTEDJOINT Joint directions, physical retraction, actual-cost acceptance.
clock=tic; candidate=seed; best=baseline; radius=c.projected_step_radius;
if nargin<6, cached=[]; end
info=struct('status','screened_feasible','exitflag',0,'failure_reason','', ...
 'max_scaled_violation',NaN,'elapsed_s',0,'replay',baseline, ...
 'selected_source','retained_feasible_seed','active_target_ids',c.joint_target_ids, ...
 'restoration',[],'fuel_started',true,'seed_visit_count',baseline.visit_count, ...
 'trials',{{}},'accepted_steps',0,'refinement_performed',false,'method','joint_sqp_physical_retraction', ...
 'fuel',struct('initial',struct('change',0)));
if ~baseline.passed, info.status='prefix_restored'; info.selected_source='retained_prefix'; end
try
 for iteration=1:c.joint_iterations
  if toc(clock)>=budget || (~isempty(c.stop_file)&&isfile(c.stop_file)), break; end
  local=c; local.shooting_reltol=min(c.shooting_reltol,1e-13);
  p=ctocscreen.v3ShootingProblem(candidate,eph,local,cached); cached=[];
  if toc(clock)>=budget || (~isempty(c.stop_file)&&isfile(c.stop_file)), break; end
  [raw,fuel]=ctocscreen.v3JointDirection(p,radius);
  info.fuel=fuel; accepted=false;
  info.refinement_performed=info.refinement_performed||fuel.refinement_performed;
  if isequal(raw,p.z0)
   radius=radius/2; continue;
  end
  % Encounter alignment retracts curved target motion as well as dynamics.
  % The free-event family preserves independently movable visits and burns.
  modes={'encounter_aligned','free_events'};
  for attempt=1:6
   alpha=[1 .3 .1]; alpha=alpha(mod(attempt-1,3)+1);
   mode=modes{ceil(attempt/3)};
   if toc(clock)>=budget, break; end
   trial=p.z0+alpha*(raw-p.z0);
   [s,connection]=p.connect(trial,max(0,budget-toc(clock)),mode);
   info.refinement_performed=info.refinement_performed||connection.passed;
   r=struct('iteration',iteration,'alpha',alpha,'radius',radius,'connection',connection, ...
    'accepted',false,'verification',[],'actual_dv_km_s',sum(vecnorm(s.delta_v_km_s,2,2)));
   if connection.passed && r.actual_dv_km_s<best.total_dv_km_s-1e-10
    % Completion of this gate is mandatory; the time limit is soft.
    check=ctocscreen.v3Verify(s,eph,c); r.verification=check;
    verified=check.passed || (numel(c.joint_target_ids)<35 && ~check.stability_repeated && ...
     ctocscreen.v3RecoveryMerit(check,c.joint_target_ids,1)==0);
    if verified&&check.visit_count>=best.visit_count
     s.witness_times_s=check.witness_times_s; s.witness_times_s(check.distance_km>1)=NaN;
     candidate=s; best=check;
     r.accepted=true; accepted=true; info.accepted_steps=info.accepted_steps+1;
     info.selected_source='solver_candidate';
    end
   end
   info.trials{end+1}=r;
   if accepted, break; end
  end
  if accepted, radius=min(c.projected_step_radius*4,radius*1.5);
  else, radius=radius/2; end
  if radius<1e-7, break; end
 end
catch err
 info.failure_reason=[err.identifier ': ' err.message];
end
info.replay=best; info.elapsed_s=toc(clock); info.budget_overrun_s=max(0,info.elapsed_s-budget);
if best.passed, info.status='screened_feasible'; else, info.status='prefix_restored'; end
% The selected replay owns the evidence, including an unchanged seed.
candidate.validation_level=best.validation_level;
end
