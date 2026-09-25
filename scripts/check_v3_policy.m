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
% --- A1: monotone, saturating plane penalty (was clipped at 0.25) ---
% Fixture note: inclination is changed by a burn along the orbit normal (here
% +z), and the magnitude is inverse-solved so the resulting TRUE inclination
% change equals the requested value -- that true change is what gets penalised.
eq=[7000 0 0 0 7.5 7.5*sind(1/2)];
vmag=norm(eq(4:6));
one=@(di)ctocscreen.v3PlaneChange(eq,[0;0;2*vmag*sind(di/2)],c).penalty;
f=@(di)arrayfun(one,di);
sample=[5.5 6 8 10 15 20 30 45 90];
assert(all(diff(f(sample))>0));                      % strictly increasing
assert(all(f(sample)<c.plane_penalty_saturation));   % bounded, never a clip jump
assert(f(10)>0.9&&f(10)<1.1);                        % recalibrated 10 deg point
assert(f(5)==0);                                     % threshold still exact
% --- A1: cost penalty monotone in the overspend, never clipped to one value ---
% Every fixture below stays under the 5 deg plane threshold so the numbers
% isolate the COST penalty. The burn has to be the one crossing the cap for the
% credit-assignment rule to charge it, hence node.J = 6.0 (< cap 6.1).
plainNode=@(J)struct('state',eq,'t',1e4,'visited',true(35,1),'J',J);
m=@()containers.Map('KeyType','char','ValueType','double');
one=@(J,dv,tbl)ctocscreen.v3PolicyFeedback(tbl,plainNode(J),0,600,dv,c);
a=one(6.0,[0 0 .5],m()); b=one(6.0,[0 0 1.2],m()); d=one(6.0,[0 0 1.3],m());
assert(a.crossing_burn&&b.crossing_burn&&d.crossing_burn);
assert(a.cost_penalty>0&&a.cost_penalty<b.cost_penalty&&b.cost_penalty<d.cost_penalty);
assert(a.penalty<b.penalty&&b.penalty<d.penalty);
assert(d.penalty<c.policy_penalty_ceiling);          % no early clipping
% Same everything except the spend already burned: reuse shares the overspend.
low=one(6.0,[0 0 1.2],m()); high=one(6.08,[0 0 1.2],m());
assert(low.cost_penalty<high.cost_penalty);
% Old code produced 4.0 for every one of these; the cost term must now differ.
% The comparison uses small plane-neutral burns because a multi-km/s burn is
% itself a large plane change and would saturate the total on its own.
% The harness derives the prefix cost from the intended full-candidate cost, so
% the rejectedCost invariant (prefix + this burn) cannot be violated by mistake.
q=@(full,dv,tbl)ctocscreen.v3PolicyFeedback(tbl, ...
 plainNode(full-norm(dv)),0,600,dv,c,full);
p20=q(6.2, [0 0 0.4],m());      % this burn crosses the cap (5.8 -> 6.2 > 6.1)
p31=q(20.9,[0 0 0.4],m());      % same burn size, cap crossed much earlier
assert(p20.crossing_burn&&~p31.crossing_burn,'A2: crossing flag wrong');
assert(p20.credit>p31.credit&&p31.credit>0,'A2: credit is not graded');
assert(p20.cost_penalty>p31.cost_penalty,'A2: later burn must not outweigh the culprit');
% A burn still inside the budget is charged nothing at all.
inside=q(5.4,[0 0 .4],m());
assert(inside.cost_penalty==0&&inside.credit==0,'A2: inside-budget burn was charged');
% Monotone in the overspend. Both fixtures keep the burn as the crossing control,
% so credit is 1 for both and only the ending spend differs.
small=q(6.3,[0 0 0.5],m());      % 5.8 -> 6.3, excess 0.2
large=q(7.0,[0 0 1.2],m());      % 5.8 -> 7.0, excess 0.9, burn still crosses
assert(small.crossing_burn&&large.crossing_burn,'A1: fixture does not cross the cap');
assert(small.credit==1&&large.credit==1,'A1: fixture credit is not 1');
assert(small.cost_penalty<large.cost_penalty&&small.penalty<large.penalty,'A1: penalty not monotone in overspend');
% Honest known property: a severely overspent candidate whose own culprit burn IS
% the large plane change saturates at the ceiling. Do not hide it.
sat=q(31,[0 0 6],m());
assert(sat.penalty<=c.policy_penalty_ceiling,'A1: penalty escaped its ceiling');
assert(sat.penalty>p31.penalty,'A1: severe candidate must rank above a mild one');
% Consistency guard rejects a prefix cost passed as the whole-candidate cost.
guardFailed=false;
try
 ctocscreen.v3PolicyFeedback(m(),plainNode(6.0),0,600,[0 0 1.2],c,6.0);
catch err
 guardFailed=strcmp(err.identifier,'ctocscreen:v3:policyFeedback');
end
assert(guardFailed,'A1: prefix-cost guard did not fire');
% --- A2: graded attribution ---
% The crossing burn carries full credit. A burn in a candidate whose cap was
% already crossed earlier still carries its own small share (so mid-course waste
% stays learnable) but strictly less credit than the crossing burn.
culprit=ctocscreen.v3PolicyFeedback(m(),plainNode(6.0),0,600,[0 0 .5],c,7.2);
bystander=ctocscreen.v3PolicyFeedback(m(),plainNode(7.5),0,600,[0 0 .5],c,8.0);
assert(culprit.crossing_burn&&~bystander.crossing_burn,'A2: crossover mislabelled');
assert(culprit.credit>0&&bystander.credit>0,'A2: credit switched off');
assert(culprit.credit>bystander.credit,'A2: late burn must carry weaker credit');
assert(culprit.cost_penalty>0&&bystander.cost_penalty>0,'A2: charge missing');
% An under-cap candidate is the only one charged nothing: here the candidate
% ends at 5.8 km/s, so the burn is not crossing the cap.
underCap=ctocscreen.v3PolicyFeedback(m(),plainNode(5.3),0,600,[0 0 .5],c,5.8);
assert(underCap.credit==0&&underCap.cost_penalty==0,'A2: under-cap was charged');
% --- A2: same action, same state, same table => no double punishment ---
replayTable=m();
first=ctocscreen.v3PolicyFeedback(replayTable,plainNode(6.0),0,600,[0 0 .5],c,6.5);
again=ctocscreen.v3PolicyFeedback(replayTable,plainNode(6.0),0,600,[0 0 .5],c,6.5);
assert(first.new_penalty>0&&again.new_penalty==0,'A2: duplicate observation was charged twice');
fprintf('POLICY_CHECK: actual inclination, complete/partial cost semantics, negative weight, free-action key, route credit, terminal/truncated trace, deduplication, decay, operator rotation, monotone saturating plane penalty, monotone overspend penalty, crossing-burn credit assignment passed.\n');
end
