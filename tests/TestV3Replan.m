classdef TestV3Replan < matlab.unittest.TestCase
 properties
  Eph
  Config
  Schedule
 end
 methods(TestClassSetup)
  function fixture(test)
   [test.Eph,test.Config,test.Schedule]=v3SyntheticFixture('separated_visits');
   test.Config.action_durations_s=[100 300]; test.Config.completion_steps=2;
  end
 end
 methods(Test)
  function probeSeparatesPhysicalWitnessFromFailedPlan(test)
   s=test.Schedule; s.visit_plan_times_s=600*ones(35,1);
   probe=ctocscreen.v3ReplanProbe(s,test.Eph,test.Config,2);
   test.verifyEqual(probe.status,'complete');
   test.verifyLessThan(max(probe.distance_km),.01);
   test.verifyGreaterThan(min(probe.planning_distance_km),1);
   [~,report]=ctocscreen.v3Replan(s,test.Eph,test.Config,stream(),tableMap(),2);
   test.verifyEqual(report.plan_miss_targets,1:35);
   test.verifyEmpty(report.unconfirmed_targets);
   test.verifyFalse(report.full_feasibility_verified);
  end
  function cachedPrefixRetainsIncomingState(test)
   s=test.Schedule; s.maneuver_times_s=[100;250]; s.delta_v_km_s=[.03 0 0;0 .04 0];
   probe=ctocscreen.v3ReplanProbe(s,test.Eph,test.Config,2);
   n=ctocscreen.v3PrefixNode(s,250,test.Eph,test.Config,probe);
   [~,trace]=ctocscreen.v3Replay(s,test.Eph,test.Config,false); x=deval(trace.arcs{2},250);
   test.verifyEqual(n.state,x(1:6).','AbsTol',1e-8);
   test.verifyEqual(n.schedule.maneuver_times_s,100,'AbsTol',1e-12);
   test.verifyEqual(n.schedule.delta_v_km_s,[.03 0 0],'AbsTol',1e-12);
  end
  function rejectsProbeFromDifferentPhysicalSchedule(test)
   probe=ctocscreen.v3ReplanProbe(test.Schedule,test.Eph,test.Config,2);
   s=test.Schedule; s.initial_q(6)=s.initial_q(6)+.001;
   test.verifyError(@()ctocscreen.v3PrefixNode(s,200,test.Eph,test.Config,probe), ...
    'ctocscreen:v3:replanProbe');
  end
  function exhaustedBudgetStartsNoPropagation(test)
   probe=ctocscreen.v3ReplanProbe(test.Schedule,test.Eph,test.Config,0);
   test.verifyEqual(probe.status,'budget_exhausted'); test.verifyEmpty(probe.trace.arcs);
   [p,report]=ctocscreen.v3Replan(test.Schedule,test.Eph,test.Config,stream(),tableMap(),0);
   test.verifyEmpty(p); test.verifyEmpty(report.attempts);
   test.verifyEqual(report.reason,'budget_exhausted');
  end
  function repeatedAttemptsRotateExpensiveCuts(test)
   s=test.Schedule; s.maneuver_times_s=[590;600]; s.delta_v_km_s=[.2 0 0;.1 0 0];
   c=test.Config; c.replan_cut_count=2;
   [~,first]=ctocscreen.v3Replan(s,test.Eph,c,stream(),tableMap(),2,1);
   [~,second]=ctocscreen.v3Replan(s,test.Eph,c,stream(),tableMap(),2,2);
   test.verifyEqual(first.cut_times_s,[590 600],'AbsTol',1e-12);
   test.verifyEqual(second.cut_times_s,[600 0],'AbsTol',1e-12);
   test.verifyEqual(second.attempt_number,2);
  end
  function canRemoveBurnAtMissionHorizon(test)
   s=test.Schedule; s.maneuver_times_s=600; s.delta_v_km_s=[.1 0 0];
   c=test.Config; c.replan_cut_count=1;
   [p,report]=ctocscreen.v3Replan(s,test.Eph,c,stream(),tableMap(),2);
   test.assertNotEmpty(p); test.verifyEqual(report.cut_times_s,600,'AbsTol',1e-12);
   test.verifyEmpty(p{1}.schedule.maneuver_times_s);
   result=ctocscreen.v3Verify(p{1}.schedule,test.Eph,c);
   test.verifyTrue(result.passed); test.verifyEqual(result.total_dv_km_s,0,'AbsTol',1e-12);
   test.verifyFalse(report.full_feasibility_verified);
  end
  function rebuiltPrefixPreservesOriginalSuffixControls(test)
   s=test.Schedule; s.maneuver_times_s=[400;550]; s.delta_v_km_s=[.001 0 0;0 .002 0];
   s.initial_q(6)=s.initial_q(6)+.0005; s.witness_times_s(:)=NaN;
   c=test.Config; c.completion_steps=1; c.action_durations_s=100; c.replan_cut_count=1;
   [p,~]=ctocscreen.v3Replan(s,test.Eph,c,stream(),tableMap(),3);
   short=p(cellfun(@(v)v.schedule.replan_origin.rebuilt_end_s<400,p)); test.assertNotEmpty(short);
   candidate=short{1}.schedule; suffix=candidate.maneuver_times_s>=400;
   test.verifyEqual(candidate.maneuver_times_s(suffix),s.maneuver_times_s,'AbsTol',1e-12);
   test.verifyEqual(candidate.delta_v_km_s(suffix,:),s.delta_v_km_s,'AbsTol',1e-12);
   test.verifyTrue(all(isnan(candidate.witness_times_s)| ...
    candidate.witness_times_s<=candidate.replan_origin.rebuilt_end_s));
   test.verifyEqual(short{1}.target_ids,1:35);
   test.verifyEqual(candidate.validation_level,'proposal');
  end
  function unscreenedCompletionCannotClaimPrefixVerification(test)
   n=ctocscreen.v3PrefixNode(test.Schedule,0,test.Eph,test.Config);
   c=test.Config; c.completion_seconds=1; c.completion_steps=1;
   [~,info,~,nodes]=ctocscreen.v3Complete(n,test.Eph,c,stream(),tableMap(),false);
   test.assertNotEmpty(nodes); test.verifyFalse(info.prefix_screening);
   test.verifyEqual(info.consistency_checks,0);
   test.verifyTrue(all(cellfun(@(n)~n.completion_check.passed,nodes)));
   test.verifyEqual(info.replay_status,'proposal_pending_full_recovery');
  end
  function failedReplanIsRetriedWithoutLosingOriginalParent(test)
   [state,parent,key]=retryFixture(test.Eph,test.Config,test.Schedule);
   attempts=state.replan_history(cellfun(@(h)strcmp(h.input_key,key),state.replan_history));
   test.verifyEqual(numel(attempts),state.config.replan_retry_limit);
   test.verifyEqual(cellfun(@(h)h.attempt_number,attempts),1:state.config.replan_retry_limit);
   test.verifyTrue(any(strcmp(state.replanned_keys,ctocscreen.v3RecoveryKey(parent.schedule,1:35,parent))));
   test.verifyEqual(cellfun(@(h)h.parent_J_km_s,attempts),repmat(parent.J,1,numel(attempts)),'AbsTol',1e-12);
   test.verifyTrue(all(cellfun(@(h)strcmp(h.parent_key,ctocscreen.v3ScheduleKey(parent.schedule)),attempts)));
   test.verifyTrue(all(cellfun(@(h)isequal(h.target_ids,1:35),attempts)));
   test.verifyTrue(all(cellfun(@(h)isempty(h.attempts),attempts)));
  end
 end
end

function r=stream()
r=RandStream('mt19937ar','Seed',1);
end
function t=tableMap()
t=containers.Map('KeyType','char','ValueType','double');
end
function [state,parent,key]=retryFixture(eph,c,s)
c.root_count=1; c.colony_count=1; c.budget_s=8; c.joint_seconds=.001;
c.completion_seconds=.001; c.beam_seconds=.001; c.completion_steps=1;
c.lineage_stall_limit=1000; c.replan_seconds=realmin; c.replan_retry_limit=2;
c.stop_file=[mfilename('fullpath') '.m'];
sim=fileparts(fileparts(mfilename('fullpath'))); base=fullfile(sim,'runs','v3','development');
folder=tempname(base); state=ctocscreen.v3Search(folder,eph,c);
parent=struct('schedule',s,'J',10); key=ctocscreen.v3ScheduleKey(s);
state.replan_queue={struct('schedule',s,'target_ids',1:35,'parent',parent,'stalls',c.recovery_stall_limit)};
c.resume_file=fullfile(folder,'resume_input.mat'); save(c.resume_file,'state');
c.stop_file=''; state=ctocscreen.v3Search(tempname(base),eph,c);
end
