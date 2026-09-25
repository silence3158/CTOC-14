function check_v3_policy()
%CHECK_V3_POLICY Focused invariants for user policy and operator scheduling.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
c=ctocscreen.v3Defaults(); table=containers.Map('KeyType','char','ValueType','double');
visited=true(35,1); visited(1)=false;
node=struct('state',[7000 0 0 0 7.5 0],'t',0,'visited',visited,'J',6);
dv=[0 0 1]; p=ctocscreen.v3PlaneChange(node.state,dv,c);
assert(p.inclination_change_deg>5&&p.penalty>0);
r=ctocscreen.v3PolicyFeedback(table,node,1,1000,dv,c);
assert(r.rejected_cost&&r.new_penalty>0);
[v,~,detail]=ctocscreen.v3Pheromone('get',table,r.policy_key,0,c);
assert(v<c.pheromone_baseline&&v>=c.pheromone_floor&&detail.negative_penalty>0);
again=ctocscreen.v3PolicyFeedback(table,node,1,1000,dv,c); assert(again.new_penalty==0);
ctocscreen.v3Pheromone('evaporate',table,'',0,c);
after=ctocscreen.v3Pheromone('get',table,r.policy_key,0,c); assert(after>v);
partial=node; partial.visited(:)=false;
coplanar=ctocscreen.v3PolicyFeedback(table,partial,2,2000,[0 1 0],c);
assert(~coplanar.rejected_cost&&coplanar.penalty==0);
free=ctocscreen.v3PolicyFeedback(table,partial,0,1000,dv,c);
assert(strcmp(free.policy_key,ctocscreen.v3PolicyKey(partial,2,0,1000))&&free.new_penalty>0);
route=ctocscreen.v3PolicyFeedback(table,partial,2,2000,[0 1 0],c,10);
assert(route.rejected_cost&&route.new_penalty>0);
rejected=ctocscreen.v3PolicyFeedback(table,partial,3,2000,[0 .5 0],c,10,false);
assert(rejected.rejected_prefix&&rejected.rejected_cost&&~rejected.complete_candidate&&rejected.new_penalty>0);
trace=struct('times',[0;1000],'states',[7000 0 0 0 7.5 0;7000 0 0 0 7.5 1]);
schedule=struct('maneuver_times_s',1000,'delta_v_km_s',[0 0 1]);
assert(abs(ctocscreen.v3TracePlanePenalty(schedule,trace,c)-p.penalty)<1e-12);
schedule.maneuver_times_s=zeros(0,1); schedule.delta_v_km_s=zeros(0,3);
assert(ctocscreen.v3TracePlanePenalty(schedule,trace,c)==0);
s=struct('root_key','same_lineage'); state=struct('work_lineage_keys',{{}},'work_lineage_attempts',[]);
ops=cell(1,6);
for k=1:6, [state,ops{k}]=ctocscreen.v3WorkOperator(state,s,c); end
assert(isequal(ops,{'joint_refine','directed_replan','mutation','joint_refine','directed_replan','mutation'}));
fprintf('POLICY_CHECK: actual inclination, complete/partial cost semantics, negative weight, free-action key, route credit, terminal/truncated trace, deduplication, decay, operator rotation passed.\n');
end
