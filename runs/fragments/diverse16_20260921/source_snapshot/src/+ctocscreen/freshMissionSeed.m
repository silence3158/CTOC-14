function [source,log]=freshMissionSeed(base,p,cfg,seed,mode)
%FRESHMISSIONSEED Construct a whole-mission seed independently of local SQP noise.
stream=RandStream('mt19937ar','Seed',seed);source=[];
log=struct('seed',seed,'mode',mode,'attempts',0,'failures',{{}});start=tic;
if strcmp(mode,'greedy')
 g=ctocscreen.defaultConfig();g.geometry_seed_probability=1;g.construct_refine_count=1;
 g.construct_time_budget_factor=1.8;g.construct_times_s=[3600 7200 14400 28800];
 g.construction_min_altitude_after_first_km=200;log.construction_config=g;
 for attempt=1:3
  log.attempts=attempt;
  try
   [c,diagnostic]=ctocscreen.constructGreedy(p,g,stream);assert(~isempty(c),diagnostic.status);
   e=ctocscreen.evaluate(c,p,g,true);assert(strcmp(e.status,'two_body_verified'));
   s=ctocscreen.importV1Candidate(c,e);nr=ctocscreen.propagateSchedule(s,p,false);assert(nr.passed);
   ir=ctocscreen.propagateSchedule(s,p,true);assert(ir.passed);
   source=struct('schedule',s,'evaluation',nr,'independent',ir);break;
  catch err,log.failures{end+1}=err.message;
  end
 end
end
if isempty(source)
 cc=cfg;cc.kick_attempts=50;cc.kick_seconds=20;
 [source,klog]=ctocscreen.kickArcSeed(base,p,cc,stream);log.kick=klog;
 if strcmp(mode,'greedy'),log.fallback='nonlocal_kick_after_failed_construction';end
end
assert(~isempty(source),'No independently verified seed was produced.');
log.elapsed_s=toc(start);log.random_state=stream.State;
log.order_hamming=sum(source.schedule.event_target_ids~=base.schedule.event_target_ids);
log.time_rms_s=sqrt(mean((source.schedule.event_times_s-base.schedule.event_times_s).^2));
log.initial_q=source.schedule.initial_q;log.initial_cost=source.independent.total_dv_km_s;
end
