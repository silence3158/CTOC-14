function result=analyze_v3_tail27(initialSource)
%ANALYZE_V3_TAIL27 Diagnostic reconstruction, never a cold-search input.
if nargin<1, initialSource='configured_start'; end
assert(ismember(initialSource,{'configured_start','cached_coast'}));
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
label='tail27_20260925';
if strcmp(initialSource,'cached_coast'), label=[label '_cached_coast']; end
folder=fullfile(sim,'runs/v3/diagnostics',label);
if ~isfolder(folder), mkdir(folder); end
assert(~isfile(fullfile(folder,'comparison.mat')),'Do not overwrite diagnostic evidence.');
doc=xmlread(fullfile(sim,'tail27_polished.atk'));
xp=javax.xml.xpath.XPathFactory.newInstance().newXPath();
seq=one(doc,'//Satellite[@Name="Inspector"]/Orbit/Propagator/Planning/Segment');
nodes=xp.evaluate('Segment',seq,javax.xml.xpath.XPathConstants.NODESET);
start=one(seq,'Segment[@ComponentType="CMCSInitialState"]/InitialState');
if strcmp(initialSource,'cached_coast')
 start=one(seq,'Segment[@Name="Coast1"]/InitialState');
end
x0=state(start); eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim, ...
 'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v3Defaults(); c.guided_max_revolutions=0; c.guided_branches=2;
c.plane_penalty_km_s=0; c.budget_s=30;
plan=zeros(35,5); time=0; coast=0; k=0; cachedDv=zeros(35,1);
for j=0:nodes.getLength()-1
 node=nodes.item(j); type=char(node.getAttribute('ComponentType'));
 assert(strcmp(txt(node,'Active'),'1'),'Inactive segment needs explicit handling.');
 switch type
  case 'CMCSInitialState'
   assert(j==0,'Only one initial state at sequence start is supported.');
  case 'CMCSPropagate'
   dt=str2double(attr(node,'StopCondition/StateTrigger','TripValue'));
   coast=coast+dt; time=time+dt;
  case 'CMCSLambertTarget'
   assert(strcmp(txt(node,'IsOnedv'),'1')&&strcmp(txt(node,'IsPerturb'),'1'));
   assert(strcmp(txt(node,'Revolution'),'0'),'This diagnostic expects zero-revolution legs.');
   k=k+1; id=sscanf(attr(node,'CoordSystem/CoordAxes/CoordPoint','Name'),'Target%d');
   dt=str2double(txt(node,'Duration')); plan(k,:)=[id,coast,time,dt,time+dt];
   cachedDv(k)=str2double(attr(node,'DV1','DV')); coast=0; time=time+dt;
  otherwise
   error('Unsupported segment type %s',type);
 end
end
assert(k==35&&numel(unique(plan(:,1)))==35);
targetDelta=zeros(35,1);
for j=1:35
 o=one(doc,sprintf('//Satellite[@Name="Target%d"]/Orbit',j));
 vals=zeros(1,6); names={'PositionX','PositionY','PositionZ','VelocityX','VelocityY','VelocityZ'};
 for a=1:6, vals(a)=str2double(txt(o,names{a}))/1000; end
 targetDelta(j)=max(abs(vals-eph.states0(j,:)));
end
result=struct('purpose','diagnostic_only_not_cold_search_or_original_pulse_verification', ...
 'initial_source',initialSource,'atk_sha256',ctocscreen.v3FileHash(fullfile(sim,'tail27_polished.atk')), ...
 'branch_rule','zero-revolution minimum departure cost; ATK branch equivalence unverified', ...
 'plan',plan,'cached_dv',cachedDv,'target_state_max_component_difference',max(targetDelta), ...
 'initial_state',x0,'elapsed_s',0);
s=struct('schema_version','free_maneuver_v3','dynamics_id','central_j2', ...
 'initial_q',toQ(x0,eph.model.mu),'duration_s',time,'maneuver_times_s',plan(:,3), ...
 'delta_v_km_s',zeros(35,3),'witness_times_s',nan(35,1));
s.witness_times_s(plan(:,1))=plan(:,5); x=x0; t=0; stats=zeros(35,11); clock=tic;
fprintf('PLAN duration %.6f d coast %.6f d; initial q %s; target difference %.3g\n', ...
 time/86400,sum(plan(:,2))/86400,mat2str(s.initial_q,12),max(targetDelta));
assert(norm(ctocscreen.initialState(s.initial_q,eph.model.mu,eph.model.re).'-x0)<1e-8);
for j=1:35
 ta=plan(j,3); tb=plan(j,5); id=plan(j,1);
 x=ctocscreen.v3Arc(x,t,ta,eph.model,c,false,true);
 goal=ctocscreen.v3QueryTargets(eph,id,tb,'pairs');
 [dv,detail,alternatives]=ctocscreen.v3GuidedTransfer(x,ta,tb,goal,eph.model,c,tic);
 s.delta_v_km_s(j,:)=dv.';
 change=ctocscreen.v3PlaneChange(x,dv,c);
 stats(j,:)=[j,id,norm(dv),norm(x(4:6)),norm(x(1:3)), ...
  change.inclination_change_deg,change.plane_rotation_deg,shape(x,eph.model.mu), ...
  shape([x(1:3);x(4:6)+dv],eph.model.mu)];
 x=alternatives{1}.final_state.'; t=tb;
 fprintf('LEG %02d target %02d dv %.6f cumulative %.6f speed %.3f di %.3f plane %.3f\n', ...
 j,id,norm(dv),sum(vecnorm(s.delta_v_km_s,2,2)),stats(j,4),stats(j,6:7));
 assert(detail.successful_branches>=1);
end
result.reconstructed_schedule=s; result.reconstructed_burns=stats;
[result.reconstructed_verification,result.reconstructed_trace]=ctocscreen.v3Replay(s,eph,c,true);
fprintf('RECONSTRUCTED verified %d visits %d dv %.9f maxdistance %.6f height %.6f\n', ...
 result.reconstructed_verification.passed,result.reconstructed_verification.visit_count, ...
 result.reconstructed_verification.total_dv_km_s,max(result.reconstructed_verification.distance_km), ...
 result.reconstructed_verification.min_altitude_lower_km);
labels={'v3_policy_cold_dev_300_seed888_09','v3_beam1_cs30_cold_dev_300_seed888_01'};
result.baselines=cell(size(labels));
for j=1:numel(labels)
 old=load(fullfile(sim,'runs/v3/search',labels{j},'rejected_complete.mat'),'rejected');
 b=old.rejected; schedule=b.schedule;
 % State diagnostics only: archived independent verification stays separate.
 [screen,trace]=ctocscreen.v3Replay(schedule,eph,c,false,'prefix_witnesses');
 burns=zeros(numel(schedule.maneuver_times_s),11);
 for a=1:size(burns,1)
  ix=find(trace.times==schedule.maneuver_times_s(a),1); dv=schedule.delta_v_km_s(a,:).';
  before=trace.states(ix,:).'; before(4:6)=before(4:6)-dv;
  change=ctocscreen.v3PlaneChange(before,dv,c);
  burns(a,:)=[a,0,norm(dv),norm(before(4:6)),norm(before(1:3)), ...
   change.inclination_change_deg,change.plane_rotation_deg,shape(before,eph.model.mu), ...
   shape([before(1:3);before(4:6)+dv],eph.model.mu)];
 end
 result.baselines{j}=struct('label',labels{j},'schedule',schedule, ...
  'archived_verification',b.verification,'diagnostic_screen',screen,'burns',burns);
end
result.elapsed_s=toc(clock); result.config=c;
result.implementation_signature=ctocscreen.v3ImplementationSignature(eph);
save(fullfile(folder,'comparison.mat'),'result','-v7.3');
writematrix(plan,fullfile(folder,'teammate_plan.csv'));
writematrix(stats,fullfile(folder,'teammate_reconstructed_burns.csv'));
for j=1:numel(result.baselines)
 writematrix(result.baselines{j}.burns,fullfile(folder,sprintf('baseline_%d_burns.csv',j)));
end
fprintf('SAVED %s elapsed %.3f s\n',folder,result.elapsed_s);

 function found=one(context,path)
  found=xp.evaluate(path,context,javax.xml.xpath.XPathConstants.NODE);
  assert(~isempty(found),'XML node missing: %s',path);
 end
 function value=txt(context,path)
  selectedText=one(context,path); value=char(selectedText.getTextContent());
 end
 function value=attr(context,path,key)
  selectedAttr=one(context,path); value=char(selectedAttr.getAttribute(key));
 end
 function x=state(context)
  x=zeros(6,1); axes={'X','Y','Z'};
  for a=1:3
   x(a)=str2double(attr(context,'Pos',axes{a}))/1000;
   x(a+3)=str2double(attr(context,'Vel',axes{a}))/1000;
  end
 end
end

function q=toQ(x,mu)
r=x(1:3); v=x(4:6); h=cross(r,v); h=h/norm(h);
inc=acos(h(3)); node=cross([0;0;1],h); node=node/norm(node);
side=cross(h,node); ev=cross(v,cross(r,v))/mu-r/norm(r);
q=[1/(2/norm(r)-dot(v,v)/mu),dot(ev,node),dot(ev,side),inc, ...
 atan2(node(2),node(1)),atan2(dot(r,side),dot(r,node))];
end

function ae=shape(x,mu)
r=x(1:3); v=x(4:6); e=norm(cross(v,cross(r,v))/mu-r/norm(r));
ae=[1/(2/norm(r)-dot(v,v)/mu),e];
end
