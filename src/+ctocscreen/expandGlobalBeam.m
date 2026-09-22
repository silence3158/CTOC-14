function out=expandGlobalBeam(source,p,cfg)
%EXPANDGLOBALBEAM Search all remaining targets in depth-indexed fragment beams.
% Positive visits only; finite beam/time/branch grids remain search approximations.
stream=RandStream('mt19937ar','Seed',cfg.seed); N=size(p.states0,1);
layers=cell(1,N+1); layers{1}={}; complete={}; start=tic;
empty=struct('schema_version','fragment_schedule_v2','dynamics_id','two_body', ...
 'initial_q',[],'maneuver_times_s',zeros(0,1),'delta_v_km_s',zeros(0,3), ...
 'duration_s',0,'event_target_ids',zeros(0,1),'event_times_s',zeros(0,1));
out=struct('config',cfg,'layers',zeros(0,4),'complete',{{}},'partial',{{}}, ...
 'lambert_calls',0,'multi_trials',0,'multi_passed',0,'source_hash',ctocscreen.implementationHash());
base=source.schedule;br=ctocscreen.propagateSchedule(base,p,false);assert(br.passed);
for j=1:cfg.beam_roots
 q=source.schedule.initial_q;
 if j>1
  q(1)=p.re_km+590+20*rand(stream); ecc=.000999*sqrt(rand(stream)); arg=2*pi*rand(stream);
  q(2:3)=ecc*[cos(arg) sin(arg)];
  if j==cfg.beam_roots, q(4)=acos(2*rand(stream)-1); q(5:6)=2*pi*rand(stream,1,2);
  else, q(4)=min(pi,max(0,q(4)+.1*randn(stream)));q(5:6)=mod(q(5:6)+.15*randn(stream,1,2),2*pi);end
 end
 s=empty; s.initial_q=q;
 layers{1}{end+1}=struct('schedule',s,'state',ctocscreen.initialState(q,p.mu_km3_s2), ...
  'time',0,'mask',false(1,N),'cost',0,'last_fragment',[]);
end
for depth=0:N-1
 if toc(start)>=cfg.beam_seconds, break;end
 if cfg.beam_warm_prefixes&&depth>0
  warm=base;warm.event_target_ids=base.event_target_ids(1:depth);warm.event_times_s=base.event_times_s(1:depth);
  warm.duration_s=warm.event_times_s(end);keep=base.maneuver_times_s<warm.duration_s;
  warm.maneuver_times_s=base.maneuver_times_s(keep);warm.delta_v_km_s=base.delta_v_km_s(keep,:);
  mask=false(1,N);mask(warm.event_target_ids)=true;
  layers{depth+1}{end+1}=struct('schedule',warm,'state',br.event_states(depth,:), ...
   'time',warm.duration_s,'mask',mask,'cost',sum(vecnorm(warm.delta_v_km_s,2,2)),'last_fragment',[]);
 end
 beam=ctocscreen.selectDiverseBeam(layers{depth+1},cfg.beam_width);
 if isempty(beam), continue;end
 out.partial=beam; before=out.lambert_calls;
 for parent=1:numel(beam)
  if toc(start)>=cfg.beam_seconds, break;end
  a=beam{parent}; remaining=find(~a.mask); children={};
  nominal=(p.horizon_s-a.time)/(numel(remaining)+.5);
  for id=remaining
   if toc(start)>=cfg.beam_seconds, break;end
   options=[zeros(numel(cfg.beam_times),1),max(300,nominal*cfg.beam_times')];
   old=find(base.event_target_ids==id,1);reference=[];
   if cfg.beam_warm_prefixes&&base.event_times_s(old)>a.time+60
    b0=find(base.maneuver_times_s<base.event_times_s(old),1,'last');
    w=max(0,base.maneuver_times_s(b0)-a.time);
    options(end+1,:)=[w base.event_times_s(old)-a.time];
    reference=br.preburn_states(b0,4:6)+base.delta_v_km_s(b0,:);
   end
   for option=1:size(options,1)
    w=options(option,1);dt=options(option,2);if a.time+dt>p.horizon_s||dt-w<60,continue;end
    coast=ctocscreen.checkArc(a.state,w,p);if ~strcmp(coast.status,'ok')||coast.min_altitude_km<200,continue;end
    depart=ctocscreen.propagateTwoBody(a.state,w,p.mu_km3_s2);
    tar=ctocscreen.targetStates(p,id,a.time+dt);
    bs=ctocscreen.enumerateBranches(depart(1:3),tar(1:3),dt-w,p.mu_km3_s2,cfg.branch_policy);
    out.lambert_calls=out.lambert_calls+1;
    if isempty(bs), continue;end
    [~,ii]=sort(arrayfun(@(b)norm(b.v_depart-depart(4:6)),bs));
    chosen=ii(1:min(2,numel(ii)));
    if ~isempty(reference),[~,ref]=min(arrayfun(@(b)norm(b.v_depart-reference),bs));chosen=unique([chosen(:);ref]);end
    for b=chosen(:)'
     impulse=bs(b).v_depart-depart(4:6);
     xp=[depart(1:3) depart(4:6)+impulse]; arc=ctocscreen.checkArc(xp,dt-w,p);
     if ~strcmp(arc.status,'ok')||arc.min_altitude_km<200, continue;end
     x=ctocscreen.propagateTwoBody(xp,dt-w,p.mu_km3_s2);
     if norm(x(1:3)-tar(1:3))>.05, continue;end
     child=a; child.state=x;child.time=a.time+dt;child.mask(id)=true;
     child.schedule.maneuver_times_s(end+1,1)=a.time+w;
     child.schedule.delta_v_km_s(end+1,:)=impulse;
     child.schedule.event_target_ids(end+1,1)=id;child.schedule.event_times_s(end+1,1)=child.time;
     child.schedule.duration_s=child.time; child.cost=sum(vecnorm(child.schedule.delta_v_km_s,2,2));
     child.last_fragment=[]; children{end+1}=child; %#ok<AGROW>
    end
   end
  end
  children=ctocscreen.selectDiverseBeam(children,cfg.beam_keep_per_parent);
  layers{depth+2}=[layers{depth+2} children];
  % Trial from a good single-target arc: discover other targets on its extension.
  if ~isempty(children)&&numel(remaining)>=2&&toc(start)<cfg.beam_seconds
   seed=children{1}; first=seed.schedule.event_target_ids(end);
   burnWait=seed.schedule.maneuver_times_s(end)-a.time;
   dA=seed.time-a.time-burnWait; u=seed.schedule.delta_v_km_s(end,:);
   horizon=min(p.horizon_s,seed.time+nominal*1.7);
   ids=setdiff(remaining,first);
   near=ctocscreen.scanEncounters(seed.state,seed.time,horizon-seed.time,p,ids,600);
   near=near(near(:,3)<=cfg.candidate_gate_km & near(:,2)>seed.time+60,:);
   if ~isempty(near)
    pair=[first near(1,1)]; z=[burnWait u dA near(1,2)-seed.time];
    lim=[0 u-2 60 60;max(dA*.5,burnWait+300) u+2 dA*1.8 nominal*2;horizon zeros(1,5)];
    cc=cfg;cc.fragment_seconds=min(2,max(.01,cfg.beam_seconds-toc(start)));
    ff=ctocscreen.refineFragment(a.state,a.time,pair,z,lim,p,cc);out.multi_trials=out.multi_trials+1;
    if ff.passed
     f=struct('passed',true,'t_in_s',a.time,'maneuver_times_s',a.time+ff.z(1), ...
      'delta_v_km_s',ff.z(2:4),'event_target_ids',pair(:),'event_times_s',ff.event_times_s, ...
      'state_out',ff.state_out,'t_out_s',ff.event_times_s(end),'counts',[2]);
     addMulti(f);
    end
   end
   if size(near,1)>=2&&depth+3<=N&&toc(start)<cfg.beam_seconds
    triples=sortrows(near(1:2,:),2); pair=[first triples(:,1)'];
    zz=[burnWait max(60,.15*dA) .85*dA diff([seed.time;triples(:,2)])' reshape([.1*u;.9*u]',1,[])];
    cc=cfg; cc.multi_seconds=min(3,max(.01,cfg.beam_seconds-toc(start)));
    [f,~]=ctocscreen.refineMultiFragment(a.state,a.time,pair,[0 3],zz,p,cc,horizon);
    out.multi_trials=out.multi_trials+1; if f.passed, addMulti(f);end
   end
  end
 end
 layers{depth+2}=ctocscreen.selectDiverseBeam(layers{depth+2},cfg.beam_width);
 out.layers(end+1,:)=[depth numel(beam) out.lambert_calls-before toc(start)];
 fprintf('BEAM depth%d nodes%d calls%d elapsed%.1f\n',depth,numel(beam),out.lambert_calls,toc(start));
end
if ~isempty(layers{N+1})
 candidates=ctocscreen.selectDiverseBeam(layers{N+1},cfg.beam_width);
 for j=1:numel(candidates)
  rr=ctocscreen.propagateSchedule(candidates{j}.schedule,p,false);
  if rr.passed, complete{end+1}=struct('schedule',candidates{j}.schedule,'evaluation',rr);end %#ok<AGROW>
 end
end
out.complete=complete;out.elapsed_s=toc(start);out.random_state=stream.State;
out.complete_independent=cell(size(complete));
out.best_verified=source;out.best_verified_dv=source.independent.total_dv_km_s;
for j=1:numel(complete)
 rr=ctocscreen.propagateSchedule(complete{j}.schedule,p,true);out.complete_independent{j}=rr;
 if rr.passed&&rr.total_dv_km_s<out.best_verified_dv-1e-8
  out.best_verified=complete{j};out.best_verified.independent=rr;out.best_verified_dv=rr.total_dv_km_s;
 end
end
 function addMulti(f)
  child=a; child.schedule.maneuver_times_s=[a.schedule.maneuver_times_s;f.maneuver_times_s];
  child.schedule.delta_v_km_s=[a.schedule.delta_v_km_s;f.delta_v_km_s];
  child.schedule.event_target_ids=[a.schedule.event_target_ids;f.event_target_ids];
  child.schedule.event_times_s=[a.schedule.event_times_s;f.event_times_s];
  child.schedule.duration_s=f.t_out_s;child.state=f.state_out;child.time=f.t_out_s;
  child.mask(f.event_target_ids)=true;child.cost=sum(vecnorm(child.schedule.delta_v_km_s,2,2));child.last_fragment=f;
  slot=sum(child.mask)+1;layers{slot}=[layers{slot} {child}];out.multi_passed=out.multi_passed+1;
 end
end
