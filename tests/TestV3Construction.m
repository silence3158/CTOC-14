classdef TestV3Construction < matlab.unittest.TestCase
 properties
  Eph
  Config
  Schedule
 end
 methods(TestClassSetup)
  function fixture(test)
   [test.Eph,test.Config,test.Schedule]=v3SyntheticFixture('separated_visits');
   test.Config.action_durations_s=[100 300 550]; test.Config.budget_s=6;
   test.Config.completion_seconds=3; test.Config.beam_diversity_scales=[25 .005 10];
  end
 end
 methods(Test)
  function estimatorUsesFutureInspectorAndPositionOnly(test)
   [eph,c,s]=v3SyntheticFixture(); c.action_durations_s=300;
   n=ctocscreen.v3PrefixNode(s,0,eph,c);
   [pairs,report]=ctocscreen.v3TimeCandidates(n,eph,c,stream(),tableMap(),tic);
   test.assertNotEmpty(pairs);
   test.verifyLessThan(max(cellfun(@(p)p.estimate,pairs)),1e-7);
   test.verifyEqual(report.method,'j2_stm_position_only');
   test.verifyEqual(report.pairs_evaluated,35);
  end
  function screeningEnumeratesMultipleTimes(test)
   n=ctocscreen.v3PrefixNode(test.Schedule,0,test.Eph,test.Config);
   [pairs,report]=ctocscreen.v3TimeCandidates(n,test.Eph,test.Config,stream(),tableMap(),tic);
   test.verifyGreaterThanOrEqual(report.propagations,3);
   test.verifyEqual(numel(unique(cellfun(@(p)p.target,pairs))),35);
   test.verifyGreaterThanOrEqual(numel(unique(cellfun(@(p)p.dt,pairs))),3);
   test.verifyLessThanOrEqual(max(cellfun(@(p)p.dt,pairs)),test.Eph.model.horizon_s);
  end
  function timeSearchRetainsAlternativesAndLocalTrials(test)
   n=ctocscreen.v3PrefixNode(test.Schedule,0,test.Eph,test.Config);
   [actions,report]=ctocscreen.v3TimeSearch(n,test.Eph,test.Config,stream(),tableMap(),tic);
   test.verifyGreaterThan(numel(actions),1);
   test.verifyTrue(any(cellfun(@(r)strcmp(r.stage,'local_time')&&r.passed,report.trials)));
   test.verifyGreaterThan(numel(unique(cellfun(@(a)a.dt,actions))),1);
   test.verifyLessThan(actionEndpointError(actions,n,test.Eph,test.Config),.051);
  end
  function localTimeSearchImprovesActualImpulse(test)
   c=test.Config; c.action_durations_s=300; c.guided_attempts=1;
   c.time_refine_iterations=3; c.time_refine_fraction=.4;
   n=ctocscreen.v3PrefixNode(test.Schedule,0,test.Eph,c); n.visited(:)=true; n.visited(1)=false;
   [~,report]=ctocscreen.v3TimeSearch(n,test.Eph,c,stream(),tableMap(),tic);
   coarse=report.trials(cellfun(@(r)strcmp(r.stage,'coarse')&&r.passed,report.trials));
   local=report.trials(cellfun(@(r)strcmp(r.stage,'local_time')&&r.passed,report.trials));
   test.assertNotEmpty(coarse); test.assertNotEmpty(local);
   test.verifyLessThan(min(cellfun(@(r)r.dv_km_s,local)),min(cellfun(@(r)r.dv_km_s,coarse)));
  end
  function guidedTransferReturnsMultiplePhysicalBranches(test)
   [branches,initial,goal,model,c]=branchFixture(test.Eph.model,test.Config);
   test.verifyGreaterThan(numel(branches),1);
   test.verifyLessThan(branchError(branches,initial,goal,model,c),.051);
   velocities=cell2mat(cellfun(@(a)a.final_state(4:6),branches,'UniformOutput',false).');
   test.verifyGreaterThan(max(vecnorm(velocities-velocities(1,:),2,2)),.1);
  end
  function expansionKeepsMultipleGuidedSuccessors(test)
   n=ctocscreen.v3PrefixNode(test.Schedule,0,test.Eph,test.Config);
   c=test.Config; c.actions_per_node=3;
   [children,~,report]=ctocscreen.v3Expand(n,test.Eph,c,stream(),tableMap(),tic);
   test.verifyGreaterThan(numel(children),3);
   test.verifyGreaterThan(report.generated,3);
   test.verifyLessThan(childrenReplayError(children,test.Eph,c),1e-7);
  end
  function completionMaintainsBeamAcrossLayers(test)
   n=ctocscreen.v3PrefixNode(test.Schedule,0,test.Eph,test.Config);
   c=test.Config; c.completion_steps=3; c.completion_seconds=8;
   [~,report,~,nodes]=ctocscreen.v3Complete(n,test.Eph,c,stream(),tableMap());
   test.verifyGreaterThan(numel(nodes),1);
   test.verifyGreaterThanOrEqual(numel(report.layer_widths),2);
   test.verifyTrue(all(report.layer_widths>1));
   test.verifyTrue(all(cellfun(@(n)n.completion_check.passed,nodes)));
   test.verifyLessThan(childrenReplayError(nodes,test.Eph,c),1e-7);
  end
  function completionSupportsMultipleRoots(test)
   a=ctocscreen.v3PrefixNode(test.Schedule,0,test.Eph,test.Config);
   s=test.Schedule; s.initial_q(6)=s.initial_q(6)+.02;
   b=ctocscreen.v3PrefixNode(s,0,test.Eph,test.Config);
   c=test.Config; c.completion_steps=1;
   [~,report,~,nodes]=ctocscreen.v3Complete({a,b},test.Eph,c,stream(),tableMap());
   test.verifyEqual(report.input_count,2);
   test.verifyGreaterThan(numel(nodes),1);
  end
  function prefixBeforeBurnKeepsRealIncomingVelocity(test)
   s=test.Schedule; s.maneuver_times_s=[100;250]; s.delta_v_km_s=[.03 0 0;0 .04 0];
   n=ctocscreen.v3PrefixNode(s,250,test.Eph,test.Config);
   [~,trace]=ctocscreen.v3Replay(s,test.Eph,test.Config,false);
   before=deval(trace.arcs{2},250);
   test.verifyEqual(n.state,before(1:6).','AbsTol',1e-8);
   test.verifyEqual(n.schedule.maneuver_times_s,100);
   test.verifyEqual(n.schedule.delta_v_km_s,[.03 0 0]);
  end
  function queuePreservesPendingCandidatesAndDropsDuplicates(test)
   a=queueItem(test.Schedule); changed=test.Schedule; changed.duration_s=590; b=queueItem(changed);
   [queue,added]=ctocscreen.v3QueueCandidates({a},{a,b},2);
   test.verifyEqual(added,1); test.verifyEqual(numel(queue),2);
   test.verifyEqual(queue{1}.schedule,a.schedule); test.verifyEqual(queue{2}.schedule,b.schedule);
   changed.duration_s=580;
   [full,added]=ctocscreen.v3QueueCandidates(queue,{queueItem(changed)},2);
   test.verifyEqual(added,0); test.verifyEqual(full,queue);
  end
  function budgetDeferralPreservesProposalEvenWhenQueueIsFull(test)
   state=struct('pending_initial',{{}},'fresh_queue',{{}},'recovery_pool',{{queueItem(test.Schedule)}});
   parent=struct('schedule',test.Schedule,'J',10);
   item=struct('schedule',test.Schedule,'target_ids',1:35,'parent',parent,'stalls',2);
   result=ctocscreen.v3DeferProposal(state,'work',item);
   test.verifyEqual(numel(result.recovery_pool),2);
   test.verifyEqual(result.recovery_pool{1},item);
   test.verifyEqual(result.recovery_pool{2},state.recovery_pool{1});
   result=ctocscreen.v3DeferProposal(state,'initial',item);
   test.verifyEqual(result.pending_initial,{test.Schedule});
   fresh=queueItem(test.Schedule); result=ctocscreen.v3DeferProposal(state,'fresh',fresh);
   test.verifyEqual(result.fresh_queue,{fresh});
  end
  function directedReplanIsEnabledByDefault(test)
   c=ctocscreen.v3Defaults(); test.verifyTrue(c.replan_enabled);
  end
 end
 methods(Test)
  function costlyConnectionReplanCanReachIndependentGate(test)
   s=test.Schedule; s.maneuver_times_s=600; s.delta_v_km_s=[.1 0 0];
   c=test.Config; c.action_durations_s=600; c.replan_cut_count=1;
   [proposals,report]=ctocscreen.v3Replan(s,test.Eph,c,stream(),tableMap(),4);
   test.assertNotEmpty(proposals,report.reason);
   checks=cellfun(@(p)ctocscreen.v3Verify(p.schedule,test.Eph,c),proposals,'UniformOutput',false);
   feasible=checks(cellfun(@(r)r.passed,checks));
   test.assertNotEmpty(feasible);
   test.verifyLessThan(min(cellfun(@(r)r.total_dv_km_s,feasible)),.1);
   test.verifyFalse(report.full_feasibility_verified);
   test.verifyEqual(report.attempts{1}.mode,'expensive_connection_replan');
  end
  function exhaustedRecoveryGetsSeparateReplanSlot(test)
   state=struct('recovery_pool',{{}},'replan_queue',{{}},'replanned_keys',{{}});
   parent=struct('schedule',test.Schedule,'J',10);
   state=ctocscreen.v3QueueRecovery(state,test.Schedule,1:35,parent,test.Config.recovery_stall_limit,test.Config);
   test.verifyEmpty(state.recovery_pool); test.assertEqual(numel(state.replan_queue),1);
   test.verifyEqual(state.replan_queue{1}.schedule,test.Schedule);
   test.verifyEqual(state.replan_queue{1}.target_ids,1:35);
   test.verifyEqual(state.replan_queue{1}.parent,parent);
   state.replanned_keys={ctocscreen.v3RecoveryKey(test.Schedule,1:35,parent)}; state.replan_queue={};
   state=ctocscreen.v3QueueRecovery(state,test.Schedule,1:35,parent,test.Config.recovery_stall_limit,test.Config);
   test.verifyEmpty(state.replan_queue);
  end
  function missingTargetsReplanKeepsAllObligations(test)
   s=test.Schedule; s.initial_q(6)=s.initial_q(6)+.0005; s.witness_times_s(:)=NaN;
   [proposals,report]=ctocscreen.v3Replan(s,test.Eph,test.Config,stream(),tableMap(),4);
   test.assertNotEmpty(proposals,report.reason); test.verifyNotEmpty(report.focus_targets);
   test.verifyTrue(all(cellfun(@(p)isequal(p.target_ids,1:35),proposals)));
   test.verifyTrue(all(cellfun(@(p)all(isfinite(p.schedule.visit_plan_times_s)),proposals)));
   c=test.Config; c.joint_seconds=.001;
   [~,info]=ctocscreen.v3JointOptimize(proposals{1}.schedule,test.Eph,c);
   test.verifyEqual(info.active_target_ids,1:35);
  end
  function searchRunsDirectedWorkWithOriginalParent(test)
   state=searchFixture(test.Eph,test.Config,test.Schedule);
   directed=state.history(cellfun(@(h)startsWith(h.operator,'directed_')&&strcmp(h.source_family,'work'),state.history));
   test.verifyNotEmpty(state.replan_history); test.assertNotEmpty(directed);
   test.verifyTrue(all(cellfun(@(h)isequal(h.task_target_ids,1:35),directed)));
   test.verifyTrue(all(cellfun(@(h)isfinite(h.parent_J_km_s),directed)));
   test.verifyTrue(state.elite.verification.passed);
   test.verifyEqual(state.signature.schema,'free_maneuver_search_v3_8');
  end
 end
 methods(Test)
  function searchConsumesMultipleConstructedCandidates(test)
   state=queuedSearchFixture(test.Eph,test.Config,test.Schedule);
   fresh=state.history(cellfun(@(h)h.from_candidate_queue,state.history));
   test.verifyGreaterThanOrEqual(numel(fresh),2);
   test.verifyGreaterThanOrEqual(numel(unique(cellfun(@(h)h.input_schedule_key,fresh,'UniformOutput',false))),2);
   test.verifyTrue(all(cellfun(@(h)numel(h.task_target_ids)>=2,fresh)));
  end
 end
end

function r=stream()
r=RandStream('mt19937ar','Seed',1);
end
function t=tableMap()
t=containers.Map('KeyType','char','ValueType','double');
end
function item=queueItem(s)
item=struct('schedule',s,'target_ids',1:35);
end
function err=actionEndpointError(actions,node,eph,c)
err=0;
for k=1:numel(actions)
 a=actions{k}; x=node.state(:); x(4:6)=x(4:6)+a.dv;
 y=ctocscreen.v3Arc(x,node.t,node.t+a.dt,eph.model,c);
 target=ctocscreen.v3QueryTargets(eph,a.target,node.t+a.dt);
 err=max(err,norm(y(1:3).'-target));
end
end
function err=childrenReplayError(nodes,eph,c)
err=0;
for k=1:numel(nodes)
 r=ctocscreen.v3Replay(nodes{k}.schedule,eph,c,false);
 err=max(err,norm(nodes{k}.state-r.final_state));
end
end
function [branches,x,goal,m,c]=branchFixture(m,c)
m.j2=0; m.horizon_s=20000; c.budget_s=10; c.guided_branches=20;
r=m.re+1000; x=[r;0;0;0;sqrt(m.mu/r);0]; goal=[0;r;0];
[~,~,branches]=ctocscreen.v3GuidedTransfer(x,0,11000,goal,m,c,tic);
end
function err=branchError(branches,x,goal,m,c)
err=0;
for k=1:numel(branches)
 y=ctocscreen.v3Arc([x(1:3);x(4:6)+branches{k}.dv],0,11000,m,c);
 err=max(err,norm(y(1:3)-goal));
end
end
function state=searchFixture(eph,c,s)
s.maneuver_times_s=600; s.delta_v_km_s=[.1 0 0];
c.initial_candidates={s}; c.budget_s=8; c.joint_seconds=.001;
c.completion_seconds=.001; c.beam_seconds=.001; c.completion_steps=1;
c.root_count=1; c.colony_count=1; c.action_durations_s=600;
c.replan_enabled=true; c.replan_probability=1; c.replan_seconds=1.5; c.replan_cut_count=1;
sim=fileparts(fileparts(mfilename('fullpath')));
base=fullfile(sim,'runs','v3','development');
c.stop_file=[mfilename('fullpath') '.m'];
folder=tempname(base); state=ctocscreen.v3Search(folder,eph,c);
% Directed work starts after a real refinement, never after a zero-step call.
joint=c; joint.joint_seconds=4; joint.joint_iterations=1; joint.stop_file='';
[~,diagnostic]=ctocscreen.v3JointOptimize(s,eph,joint);
assert(diagnostic.refinement_performed,'Fixture must perform a joint refinement.');
state.refined_keys={ctocscreen.v3ScheduleKey(s)}; state.pending_initial={};
c.resume_file=fullfile(folder,'resume_input.mat'); save(c.resume_file,'state','diagnostic');
c.stop_file=''; state=ctocscreen.v3Search(tempname(base),eph,c);
end
function state=queuedSearchFixture(eph,c,s)
c.root_count=1; c.colony_count=1; c.beam_width=4; c.completion_steps=1;
c.budget_s=6; c.completion_seconds=.5; c.beam_seconds=.01; c.beam_expansions=1;
c.joint_seconds=.001; c.prefix_joint_every=1; c.lineage_stall_limit=1000;
c.stop_file=[mfilename('fullpath') '.m'];
sim=fileparts(fileparts(mfilename('fullpath'))); base=fullfile(sim,'runs','v3','development');
folder=tempname(base); state=ctocscreen.v3Search(folder,eph,c);
node=ctocscreen.v3PrefixNode(s,0,eph,c);
state.beam={node}; state.colonies={{node}}; state.pending_roots={node};
c.resume_file=fullfile(folder,'resume_input.mat'); save(c.resume_file,'state');
c.stop_file=''; state=ctocscreen.v3Search(tempname(base),eph,c);
end
