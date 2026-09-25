classdef TestV3RecoveryFlow < matlab.unittest.TestCase
 properties
  Eph
  Config
  Schedule
 end
 methods(TestClassSetup)
  function fixture(test)
   [test.Eph,test.Config,test.Schedule]=v3SyntheticFixture('separated_visits');
  end
 end
 methods(Test)
  function defaultPrefixCadenceSurvivesBusySources(test)
   [state,due]=busySources(test.Config);
   test.verifyEqual(due,[7 17]);
   fresh=state.source_history(cellfun(@(h)strcmp(h.source,'fresh'),state.source_history));
   test.verifyEqual(cellfun(@(h)h.source_attempt,fresh),1:6);
  end
  function defaultCadenceReachesJointCallWithBothPools(test)
   state=busySearch(test.Eph,test.Config,test.Schedule);
   calls=state.history(cellfun(@(h)h.from_candidate_queue,state.history));
   test.assertNotEmpty(calls);
   test.verifyEqual(calls{1}.iteration,7);
   test.verifyEqual(state.source_history{calls{1}.source_attempt_id}.source_attempt,3);
   test.verifyEqual(state.config.prefix_joint_every,3);
  end
  function queueKeepsDifferentVisitPlans(test)
   a=test.Schedule; a.visit_plan_times_s=a.witness_times_s;
   b=a; b.visit_plan_times_s([1 2])=flipud(b.visit_plan_times_s([1 2]));
   first=struct('schedule',a,'target_ids',1:35);
   second=struct('schedule',b,'target_ids',1:35);
   [queue,added]=ctocscreen.v3QueueCandidates({first},{second,first},4);
   test.verifyEqual(added,1); test.verifyEqual(numel(queue),2);
   test.verifyEqual(queue{2}.schedule.visit_plan_times_s,b.visit_plan_times_s,'AbsTol',1e-12);
   test.verifyEqual(ctocscreen.v3ScheduleKey(a),ctocscreen.v3ScheduleKey(b));
  end
  function explicitPlanIgnoresIncidentalWitnessRelabel(test)
   a=test.Schedule; a.visit_plan_times_s=a.witness_times_s;
   b=a; b.witness_times_s(:)=NaN;
   test.verifyEqual(ctocscreen.v3RecoveryKey(a,1:35),ctocscreen.v3RecoveryKey(b,35:-1:1));
   test.verifyNotEqual(ctocscreen.v3RecoveryKey(a,1:35),ctocscreen.v3RecoveryKey(test.Schedule,1:35));
  end
  function completedReplanDoesNotSuppressDifferentPlan(test)
   a=test.Schedule; a.visit_plan_times_s=a.witness_times_s;
   b=a; b.visit_plan_times_s([1 2])=flipud(b.visit_plan_times_s([1 2]));
   parent=struct('schedule',test.Schedule,'J',10);
   state=struct('recovery_pool',{{}},'replan_queue',{{}}, ...
    'replanned_keys',{{ctocscreen.v3RecoveryKey(a,1:35,parent)}});
   state=ctocscreen.v3QueueRecovery(state,a,1:35,parent,test.Config.recovery_stall_limit,test.Config);
   test.verifyEmpty(state.replan_queue);
   state=ctocscreen.v3QueueRecovery(state,b,1:35,parent,test.Config.recovery_stall_limit,test.Config);
   test.assertEqual(numel(state.replan_queue),1);
   test.verifyEqual(state.replan_queue{1}.parent,parent);
   test.verifyEqual(state.replan_queue{1}.schedule,b);
  end
  function queueKeepsDifferentSAParents(test)
   first=struct('schedule',test.Schedule,'target_ids',1:35, ...
    'parent',struct('schedule',test.Schedule,'J',10));
   second=first; second.parent.schedule.initial_q(6)=second.parent.schedule.initial_q(6)+.001;
   [queue,added]=ctocscreen.v3QueueCandidates({first},{second},4);
   test.verifyEqual(added,1); test.verifyEqual(numel(queue),2);
   test.verifyEqual(queue{2}.parent,second.parent);
  end
  function failedAttemptsReleaseOtherSources(test)
   state=struct('source_history',{{struct('source','fresh','outcome','joint_completed')}});
   [state,a]=ctocscreen.v3BeginAttempt(state,true,true,false);
   state.source_history{end}.outcome='construction_failed';
   [state,b]=ctocscreen.v3BeginAttempt(state,true,true,false);
   state.source_history{end}.outcome='joint_failed';
   [state,d]=ctocscreen.v3BeginAttempt(state,true,true,false);
   state.source_history{end}.outcome='construction_failed';
   [state,e]=ctocscreen.v3BeginAttempt(state,true,true,false);
   test.verifyEqual({a,b,d,e},{'work','recovery','work','fresh'});
   test.verifyEqual(numel(state.source_history),5);
  end
  function deleteAlwaysChangesTopology(test)
   s=test.Schedule; s.maneuver_times_s=[200;300]; s.delta_v_km_s=[.02 0 0;0 .03 0];
   changed=ctocscreen.v3Mutate(s,test.Eph,test.Config,RandStream('mt19937ar','Seed',2),2);
   test.verifyEqual(numel(changed.maneuver_times_s),1);
  end
  function smallRandomDeleteAndMergeAreReal(test)
   s=test.Schedule; s.maneuver_times_s=[200;300]; s.delta_v_km_s=[.02 0 0;0 .03 0];
   c=test.Config; c.wide_mutation_probability=0;
   counts=randomTopologyCounts(s,test.Eph,c);
   test.verifyNotEmpty(counts);
   test.verifyEqual(counts,ones(size(counts)));
  end
  function constructorIncludesInteriorVisits(test)
   s=test.Schedule; s.maneuver_times_s=0; s.delta_v_km_s=[.02 .01 0];
   s.visit_plan_times_s=s.witness_times_s;
   [candidate,report]=ctocscreen.v3FitVisitArcs(s,test.Eph,test.Config,10);
   checked=ctocscreen.v3Verify(candidate,test.Eph,test.Config);
   test.verifyTrue(report.passed,report.reason);
   test.verifyEqual(report.arcs{1}.target_ids,1:35);
   test.verifyLessThan(max(report.arcs{1}.after_km),1e-4);
   test.verifyTrue(checked.passed);
  end
  function deletedTopologyCanReachIndependentGate(test)
   s=test.Schedule; s.maneuver_times_s=[250;251]; s.delta_v_km_s=[.02 0 0;-.02 0 0];
   [candidate,name,report]=ctocscreen.v3ConnectedMutation(s,test.Eph,test.Config,RandStream('mt19937ar','Seed',5),10,2);
   checked=ctocscreen.v3Verify(candidate,test.Eph,test.Config);
   test.verifyEqual(name,'connected_delete');
   test.verifyTrue(report.passed,report.reason);
   test.verifyLessThan(numel(candidate.maneuver_times_s),2);
   test.verifyTrue(checked.passed);
   test.verifyFalse(report.full_feasibility_verified);
  end
  function reassignReachesConstructor(test)
   s=test.Schedule; s.maneuver_times_s=0; s.delta_v_km_s=[.001 0 0];
   [candidate,name,report]=ctocscreen.v3ConnectedMutation(s,test.Eph,test.Config,RandStream('mt19937ar','Seed',9),10,5);
   test.verifyEqual(name,'connected_reassign_witness');
   test.verifyTrue(report.visit_assignment_changed);
   test.verifyNotEqual(candidate.visit_plan_times_s,s.witness_times_s);
   arcs=[report.arcs{:}];
   test.verifyEqual(sort([arcs.target_ids]),1:35);
  end
  function oldVisitsDoNotProveNewAssignment(test)
   s=test.Schedule; s.visit_plan_times_s=s.witness_times_s;
   s.visit_plan_times_s([1 35])=s.visit_plan_times_s([35 1]);
   r=ctocscreen.v3JointReplay(s,test.Eph,test.Config,true);
   test.verifyTrue(r.physical_passed);
   test.verifyFalse(r.passed);
   test.verifyGreaterThan(max(r.planning_distance_km),1);
  end
  function labelOnlyMutationIsNotPhysicalNovelty(test)
   s=test.Schedule; s.witness_times_s=flipud(s.witness_times_s);
   difference=ctocscreen.v3MutationDifference(test.Schedule,s,test.Config,test.Eph.model);
   test.verifyFalse(difference.novel);
  end
  function infeasibleRecoveryUsesMissNotFuel(test)
   a=struct('height_passed',true,'distance_km',2*ones(35,1),'total_dv_km_s',1);
   b=a; b.distance_km=1.5*ones(35,1); b.total_dv_km_s=2;
   before=ctocscreen.v3RecoveryMerit(a,1:35,.99);
   after=ctocscreen.v3RecoveryMerit(b,1:35,.99);
   test.verifyLessThan(after,before);
  end
  function fullRecoveryKeepsMissingTargets(test)
   s=test.Schedule; s.initial_q(6)=s.initial_q(6)+.001;
   s.witness_times_s(30:35)=NaN;
   c=test.Config; c.joint_seconds=3; c.joint_iterations=1;
   [~,info]=ctocscreen.v3JointOptimize(s,test.Eph,c);
   test.verifyEqual(info.active_target_ids,1:35);
   test.verifyTrue(isfinite(info.recovery_initial_merit));
   test.verifyLessThanOrEqual(info.recovery_final_merit,info.recovery_initial_merit);
  end
  function searchRecoveryRetainsFullScope(test)
   s=test.Schedule; s.initial_q(6)=s.initial_q(6)+.001; s.witness_times_s(30:35)=NaN;
   c=test.Config; c.budget_s=6; c.joint_seconds=.001; c.initial_candidates={s};
   c.cold_completion_enabled=false; % Isolate the existing full-scope NLP retry contract.
   c.root_count=1; c.colony_count=2; c.actions_per_node=1; c.completion_steps=1;
   c.completion_seconds=.001; c.beam_seconds=.001; c.action_durations_s=60;
   sim=fileparts(fileparts(mfilename('fullpath')));
   folder=tempname(fullfile(sim,'runs','v3','development'));
   state=ctocscreen.v3Search(folder,test.Eph,c);
   recovered=state.history(cellfun(@(h)strcmp(h.source_family,'recovery'),state.history));
   test.assertNotEmpty(recovered);
   test.verifyTrue(all(cellfun(@(h)isequal(h.diagnostic.active_target_ids,1:35),recovered)));
   retries=recovered(cellfun(@(h)strcmp(h.operator,'resume_recovery'),recovered));
   test.verifyNotEmpty(retries);
   test.verifyEqual(numel(state.source_history),numel(state.construction_history));
   records=[state.construction_history{:}];
   test.verifyTrue(any([records.colony]==1 & ~[records.protected_root]));
   test.verifyTrue(any([records.colony]==2 & ~[records.protected_root]));
  end
 end
end

function [state,due]=busySources(c)
state=struct('source_history',{{}}); due=[];
for k=1:18
 [state,source]=ctocscreen.v3BeginAttempt(state,true,true,false);
 state.source_history{end}.outcome='failed_or_deferred';
 if strcmp(source,'fresh')&&mod(state.source_history{end}.source_attempt,c.prefix_joint_every)==0
  due(end+1)=k;
 end
end
end

function state=busySearch(eph,c,s)
c.root_count=1; c.colony_count=1; c.initial_candidates={s};
c.budget_s=2; c.joint_seconds=.001; c.completion_seconds=.001; c.beam_seconds=.001;
c.completion_steps=1; c.replan_enabled=false; c.lineage_stall_limit=1000;
c.stop_file=[mfilename('fullpath') '.m'];
sim=fileparts(fileparts(mfilename('fullpath'))); base=fullfile(sim,'runs','v3','development');
folder=tempname(base); state=ctocscreen.v3Search(folder,eph,c);
node=ctocscreen.v3PrefixNode(s,300,eph,c);
state.fresh_queue={struct('schedule',node.schedule,'target_ids',find(node.visited).', ...
 'node',node,'operator','beam_completion')};
changed=s; changed.initial_q(6)=changed.initial_q(6)+.001;
state.recovery_pool={struct('schedule',changed,'target_ids',1:35,'parent',[],'stalls',0)};
state.pending_initial={}; state.source_history={}; state.iteration=6;
for k=1:6, [state,~]=ctocscreen.v3BeginAttempt(state,true,true,false); end
c.resume_file=fullfile(folder,'resume_input.mat'); save(c.resume_file,'state');
c.stop_file=''; state=ctocscreen.v3Search(tempname(base),eph,c);
end

function counts=randomTopologyCounts(s,eph,c)
stream=RandStream('mt19937ar','Seed',22); counts=[];
for k=1:40
 [changed,name]=ctocscreen.v3Mutate(s,eph,c,stream);
 if ismember(name,{'delete','merge'}), counts(end+1)=numel(changed.maneuver_times_s); end
end
end
