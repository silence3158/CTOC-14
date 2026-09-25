classdef TestV3Search < matlab.unittest.TestCase
 properties
  Eph
  Config
  Schedule
 end
 methods(TestClassSetup)
  function fixture(test)
   [test.Eph,test.Config,test.Schedule]=v3SyntheticFixture();
  end
 end
 methods(Test)
  function connectionUsesRealState(test)
   s=test.Schedule; s.maneuver_times_s=[100;300]; s.delta_v_km_s=[.02 0 0;0 .03 0];
   [~,trace]=ctocscreen.v3Replay(s,test.Eph,test.Config,true);
   anchors=[trace.states(2,1:3);trace.states(end,1:3)];
   % trace knot 2 is first burn; next burn is knot 3.
   anchors(1,:)=trace.states(3,1:3);
   changed=s; changed.initial_q(6)=changed.initial_q(6)+1e-4;
   [connected,report]=ctocscreen.v3ConnectControls(changed,anchors,test.Eph.model,test.Config,15);
   [~,actual]=ctocscreen.v3Replay(connected,test.Eph,test.Config,true);
   test.verifyTrue(report.passed,report.reason);
   test.verifyLessThan(norm(actual.states(3,1:3)-anchors(1,:)),2e-5);
   test.verifyLessThan(norm(actual.states(end,1:3)-anchors(2,:)),2e-5);
   test.verifyEqual(connected.initial_q,changed.initial_q);
   test.verifyEqual(connected.maneuver_times_s,changed.maneuver_times_s);
  end
  function exhaustedConnectionIsNotFeasible(test)
   s=test.Schedule; s.maneuver_times_s=100; s.delta_v_km_s=[.02 0 0];
   [~,report]=ctocscreen.v3ConnectControls(s,[7000 0 0],test.Eph.model,test.Config,0);
   test.verifyFalse(report.passed);
   test.verifySubstring(report.reason,'connectBudget');
  end
  function forceJacobian(test)
   m=test.Eph.model; r=[7100;1200;800]; t=321;
   [a,A]=ctocscreen.v3Force(t,r,m);
   d=[.3;-.7;.2]; h=.01;
   numerical=(ctocscreen.v3Force(t,r+h*d,m)-ctocscreen.v3Force(t,r-h*d,m))/(2*h);
   test.verifyLessThan(norm(numerical-A*d),1e-12);
   test.verifyEqual(a,ctocscreen.v3ReferenceForce(t,r.',m).','AbsTol',1e-15);
  end
  function stateTransition(test)
   x=ctocscreen.initialState(test.Schedule.initial_q,test.Eph.model.mu,test.Eph.model.re).';
   [~,P]=ctocscreen.v3Arc(x,0,300,test.Eph.model,test.Config,true);
   d=[.3;.2;-.1;.001;-.002;.001]; h=.001;
   yp=ctocscreen.v3Arc(x+h*d,0,300,test.Eph.model,test.Config);
   ym=ctocscreen.v3Arc(x-h*d,0,300,test.Eph.model,test.Config);
   test.verifyLessThan(norm((yp-ym)/(2*h)-P*d),2e-5);
  end
  function simultaneousImpulses(test)
   s=test.Schedule; s.maneuver_times_s=[200;200;400]; s.delta_v_km_s=[1 0 0;-1 0 0;0 .1 0];
   s=ctocscreen.v3Normalize(s,test.Eph.model);
   test.verifyEqual(s.maneuver_times_s,400); test.verifyEqual(s.delta_v_km_s,[0 .1 0]);
  end
  function coastVisitsWithoutBurns(test)
   result=ctocscreen.v3Verify(test.Schedule,test.Eph,test.Config);
   test.verifyTrue(result.passed,result.failure_reason); test.verifyEqual(result.visit_count,35);
   test.verifyEqual(result.total_dv_km_s,0); test.verifyTrue(result.height_passed);
  end
  function impulseDoesNotResetVelocity(test)
   s=test.Schedule; s.maneuver_times_s=300; s.delta_v_km_s=[.1 .2 -.1];
   [~,trace]=ctocscreen.v3Replay(s,test.Eph,test.Config,false);
   before=deval(trace.arcs{1},300);
   test.verifyEqual(trace.states(2,4:6)-before(4:6).',s.delta_v_km_s,'AbsTol',1e-12);
   test.verifyEqual(trace.states(2,1:3),before(1:3).','AbsTol',1e-12);
  end
  function altitudeViolation(test)
   m=test.Eph.model; x=[m.re+150;0;0;0;sqrt(m.mu/(m.re+150));0];
   [~,~,sol]=ctocscreen.v3Arc(x,0,60,m,test.Config);
   h=ctocscreen.v3Height(sol,m,test.Config);
   test.verifyFalse(h.passed);
  end
  function invalidInitialOrbit(test)
   s=test.Schedule; s.initial_q(1)=test.Eph.model.re+620;
   test.verifyError(@()ctocscreen.v3Normalize(s,test.Eph.model),'ctocscreen:v3:initialOrbit');
  end
  function allJointVariablesAndGradient(test)
   s=test.Schedule; s.maneuver_times_s=150; s.delta_v_km_s=[.001 .002 -.001];
   s.witness_times_s=linspace(250,450,35).';
   s.duration_s=550;
   p=ctocscreen.v3ShootingProblem(s,test.Eph,test.Config);
   [~,~,G,GE]=p.constraints(p.z0);
   direction=sin((1:numel(p.z0)).'); direction=direction/norm(direction); h=1e-6;
   [cp,ep]=p.constraints(p.z0+h*direction); [cm,em]=p.constraints(p.z0-h*direction);
   numerical=[(cp-cm);(ep-em)]/(2*h); exact=[G.';GE.']*direction;
   test.verifyLessThan(norm(numerical-exact)/max(1,norm(exact)),2e-5);
   test.verifyTrue(issparse(G)&&issparse(GE));
   test.verifyEqual(numel(p.q_indices),6);
   test.verifyGreaterThan(p.variable_count,6+4+1+35);
  end
  function zeroDurationArc(test)
   x=ctocscreen.initialState(test.Schedule.initial_q,test.Eph.model.mu,test.Eph.model.re);
   [y,P]=ctocscreen.v3Arc(x,100,100,test.Eph.model,test.Config,true);
   test.verifyEqual(y,x.'); test.verifyEqual(P,eye(6));
  end
  function multirevolutionCoast(test)
   m=test.Eph.model; m.j2=0; m.horizon_s=20000;
   x=ctocscreen.initialState(test.Schedule.initial_q,m.mu,m.re);
   y=ctocscreen.v3Arc(x,0,12000,m,test.Config);
   reference=ctocscreen.propagateTwoBody(x,12000,m.mu);
   test.verifyLessThan(norm(y(1:3)-reference(1:3).'),.001);
  end
  function hiddenInteriorLowPoint(test)
   m=test.Eph.model; m.j2=0; m.horizon_s=20000;
   ra=m.re+600; rp=m.re+100; a=(ra+rp)/2; T=2*pi*sqrt(a^3/m.mu);
   x=[ra;0;0;0;sqrt(m.mu*(2/ra-1/a));0];
   [~,~,sol]=ctocscreen.v3Arc(x,0,T,m,test.Config);
   ends=deval(sol,[0 T]); h=ctocscreen.v3Height(sol,m,test.Config);
   test.verifyGreaterThan(min(vecnorm(ends(1:3,:),2,1))-m.re,590);
   test.verifyFalse(h.passed);
  end
  function independentTargetsRejectFalseWitnesses(test)
   eph=test.Eph; eph.states0(:,1)=eph.states0(:,1)+1000;
   result=ctocscreen.v3Verify(test.Schedule,eph,test.Config);
   test.verifyFalse(result.passed); test.verifyLessThan(result.visit_count,35);
  end
  function mutationChangesEventAssignment(test)
   s=test.Schedule; s.maneuver_times_s=200; s.delta_v_km_s=[.1 0 0];
   [changed,name]=ctocscreen.v3Mutate(s,test.Eph,test.Config,RandStream('mt19937ar','Seed',9),5);
   test.verifyEqual(name,'reassign_witness');
   test.verifyNotEqual(changed.witness_times_s,s.witness_times_s);
  end
  function zeroVisitBeamQuota(test)
   roots=ctocscreen.v3Roots(test.Eph,test.Config,RandStream('mt19937ar','Seed',2));
   a=roots{1}; a.last_gain=0; a.estimate=100; b=roots{2}; b.last_gain=1; b.estimate=1;
   beam=ctocscreen.v3SelectBeam({a,b},2);
   test.verifyEqual(numel(beam),2); test.verifyTrue(any(cellfun(@(n)n.last_gain==0,beam)));
  end
  function saAllowsUphill(test)
   [accept,p,draw]=ctocscreen.v3SAAccept(6,6.01,100,RandStream('mt19937ar','Seed',1));
   test.verifyTrue(accept); test.verifyEqual(p,exp(-.01/100),'AbsTol',1e-15);
   test.verifyLessThan(draw,p);
  end
  function jointSolverRuns(test)
   [s,info]=ctocscreen.v3JointOptimize(test.Schedule,test.Eph,test.Config);
   test.verifyNotEqual(info.status,'failed',info.failure_reason);
   test.verifyTrue(isfinite(info.exitflag)); test.verifyEqual(s.schema_version,'free_maneuver_v3');
   test.verifyEqual(info.status,'screened_feasible');
  end
  function feasibleSeedSurvivesShortBudget(test)
   c=test.Config; c.joint_iterations=1;
   [~,info]=ctocscreen.v3JointOptimize(test.Schedule,test.Eph,c);
   test.verifyEqual(info.status,'screened_feasible');
   test.verifyEqual(info.selected_source,'retained_feasible_seed');
  end
  function beamEndpointIsReal(test)
   c=test.Config; c.root_count=1; c.actions_per_node=4; c.action_durations_s=300; c.budget_s=20;
   stream=RandStream('mt19937ar','Seed',15); nodes=ctocscreen.v3Roots(test.Eph,c,stream);
   n=nodes{1}; n.schedule=test.Schedule; n.state=test.Eph.states0(1,:); n.schedule.witness_times_s(:)=NaN;
   ph=containers.Map('KeyType','char','ValueType','double');
   children=ctocscreen.v3Expand(n,test.Eph,c,stream,ph,tic);
   test.assertNotEmpty(children);
   test.verifyLessThan(endpointError(children,test.Eph,c),1e-8);
  end
  function endBurnDurationMutation(test)
   s=test.Schedule; s.maneuver_times_s=s.duration_s; s.delta_v_km_s=[.01 0 0];
   s.witness_times_s(1)=NaN;
   changed=ctocscreen.v3Mutate(s,test.Eph,test.Config,RandStream('mt19937ar','Seed',21),6);
   test.verifyLessThanOrEqual(changed.maneuver_times_s,changed.duration_s);
   test.verifyTrue(isnan(changed.witness_times_s(1)));
  end
  function recoveryScheduleRemainsValid(test)
   s=test.Schedule; s.maneuver_times_s=100; s.delta_v_km_s=[.4 .3 .2];
   c=test.Config; c.joint_iterations=1;
   [changed,~]=ctocscreen.v3JointOptimize(s,test.Eph,c);
   normalized=ctocscreen.v3Normalize(changed,test.Eph.model);
   test.verifyTrue(all(normalized.maneuver_times_s<=normalized.duration_s));
  end
  function coverageSurvivesCheapCoasts(test)
   c=test.Config; c.root_count=6;
   nodes=ctocscreen.v3Roots(test.Eph,c,RandStream('mt19937ar','Seed',12));
   nodes{1}.visited(1:6)=true; nodes{1}.estimate=100; nodes{1}.last_gain=1;
   nodes{2}.estimate=1; nodes{3}.estimate=2; nodes{4}.estimate=3;
   chosen=ctocscreen.v3SelectBeam(nodes,4);
   test.verifyGreaterThanOrEqual(max(cellfun(@(n)sum(n.visited),chosen)),6);
   test.verifyTrue(any(cellfun(@(n)n.last_gain==0,chosen)));
  end
  function initialNodesAreDynamicallyContinuous(test)
   s=test.Schedule; s.initial_q(6)=s.initial_q(6)+.001;
   p=ctocscreen.v3ShootingProblem(s,test.Eph,test.Config);
   [ci,eq]=p.constraints(p.z0);
   test.verifyLessThan(max(abs(eq)),1e-8);
   test.verifyGreaterThan(max(ci),1e-4);
  end
  function restorationReducesActualMiss(test)
   s=test.Schedule; s.initial_q(6)=s.initial_q(6)+.001;
   c=test.Config; c.joint_seconds=12; c.restoration_iterations=35;
   [~,info]=ctocscreen.v3JointOptimize(s,test.Eph,c);
   test.assertNotEmpty(info.restoration,info.failure_reason);
   test.verifyLessThan(info.restoration.final_merit,info.restoration.initial_merit*1e-3);
   test.verifyEqual(info.replay.visit_count,35);
   test.verifyTrue(info.refinement_performed);
  end
  function completionOwnsItsBudget(test)
   c=test.Config; c.root_count=1; c.budget_s=1e-9; c.completion_seconds=.5;
   c.action_durations_s=300; c.actions_per_node=3;
   stream=RandStream('mt19937ar','Seed',3); nodes=ctocscreen.v3Roots(test.Eph,c,stream);
   [~,stats]=ctocscreen.v3Complete(nodes{1},test.Eph,c,stream,containers.Map('KeyType','char','ValueType','double'));
   test.verifyGreaterThanOrEqual(stats.steps,1);
  end
  function nodeBoundsPreserveInputTrajectory(test)
   c=test.Config; c.node_position_bound=100; c.node_velocity_bound=.01;
   p=ctocscreen.v3ShootingProblem(test.Schedule,test.Eph,c);
   [~,eq]=p.constraints(p.z0);
   test.verifyGreaterThan(p.expanded_bound_count,0);
   test.verifyLessThan(max(abs(eq)),1e-8);
  end
  function stablePrefixPreservesSyntheticOrbit(test)
   [s,r,report]=ctocscreen.v3StablePrefix(test.Schedule,test.Eph,test.Config);
   test.verifyTrue(report.passed); test.verifyFalse(report.truncated);
   test.verifyEqual(s.duration_s,test.Schedule.duration_s);
   test.verifyEqual(r.visit_count,35);
  end
  function rollbackRemovesOnlyMatchingExtensions(test)
   nodes=ctocscreen.v3Roots(test.Eph,test.Config,RandStream('mt19937ar','Seed',2));
   safe=nodes{1}; safe.schedule.duration_s=100;
   bad=safe; bad.schedule.duration_s=500;
   other=nodes{2}; other.schedule.duration_s=500;
   [beam,removed]=ctocscreen.v3PruneUnstableBeam({bad,other,safe},safe);
   test.verifyEqual(removed,1); test.verifyEqual(numel(beam),2);
   test.verifyEqual(beam{1}.schedule.initial_q,other.schedule.initial_q);
  end
  function freshAndRecoveryHaveSeparateSlots(test)
   selected=arrayfun(@(n)ctocscreen.v3ChooseSource(n,false,true),0:5,'UniformOutput',false);
   test.verifyEqual(selected,{'fresh','fresh','recovery','fresh','fresh','recovery'});
  end
  function missingWitnessesStayUnknownUnderMutation(test)
   s=test.Schedule; s.witness_times_s(:)=NaN;
   changed=ctocscreen.v3Mutate(s,test.Eph,test.Config,RandStream('mt19937ar','Seed',1),5);
   test.verifyTrue(all(isnan(changed.witness_times_s)));
  end
  function fuelStepHasIndependentBenefit(test)
   s=test.Schedule; s.maneuver_times_s=s.duration_s; s.delta_v_km_s=[.1 0 0];
   c=test.Config; c.joint_seconds=8; c.joint_iterations=15;
   [candidate,info]=ctocscreen.v3JointOptimize(s,test.Eph,c);
   r=ctocscreen.v3Replay(candidate,test.Eph,c,true);
   test.verifyTrue(r.passed,info.failure_reason);
   test.verifyLessThan(r.total_dv_km_s,.1);
   test.verifyEqual(info.fuel.initial.change,0);
  end
  function guidedRootRetainsFormalInitialBounds(test)
   c=test.Config; c.guided_root_fraction=1; c.root_count=1;
   nodes=ctocscreen.v3Roots(test.Eph,c,RandStream('mt19937ar','Seed',101));
   s=ctocscreen.v3Normalize(nodes{1}.schedule,test.Eph.model);
   test.verifyEqual(nodes{1}.root_method,'target_geometry');
   test.verifyGreaterThan(nodes{1}.seed_target,0);
   test.verifyLessThan(hypot(s.initial_q(2),s.initial_q(3)),.001);
   test.verifyEqual(nodes{1}.estimate,0);
  end
  function improvedPrefixReplacesStaleDescendants(test)
   nodes=ctocscreen.v3Roots(test.Eph,test.Config,RandStream('mt19937ar','Seed',4));
   old=nodes{1}; old.schedule.duration_s=100;
   descendant=old; descendant.schedule.duration_s=500;
   improved=old; improved.schedule.initial_q(6)=improved.schedule.initial_q(6)+.001;
   selected=ctocscreen.v3ReplacePrefix({old,descendant,nodes{2}},old.schedule,improved);
   test.verifyEqual(numel(selected),2);
   test.verifyEqual(selected{1}.schedule,improved.schedule);
   test.verifyEqual(selected{2}.schedule,nodes{2}.schedule);
  end
  function verifiedInitialSeedBootstrapsFullWorkPool(test)
   sim=fileparts(fileparts(mfilename('fullpath')));
   folder=tempname(fullfile(sim,'runs','v3','development'));
   c=test.Config; c.budget_s=.1; c.root_count=1; c.colony_count=1;
   c.actions_per_node=1; c.completion_steps=1; c.initial_candidates={test.Schedule};
   state=ctocscreen.v3Search(folder,test.Eph,c);
   test.assertNotEmpty(state.elite); test.verifyTrue(state.elite.verification.passed);
   test.verifyNotEmpty(state.work_pool);
   test.verifyNotEmpty(state.pheromone_keys);
   test.verifyEqual(state.dataset_kind,'synthetic_test_only');
  end
 end
end

function worst=endpointError(children,eph,c)
worst=0;
for k=1:numel(children)
 out=ctocscreen.v3Replay(children{k}.schedule,eph,c,false);
 worst=max(worst,norm(out.final_state-children{k}.state));
end
end
