classdef TestV3Diversity < matlab.unittest.TestCase
 properties
  Eph
  Config
  Schedule
  Features
 end
 methods(TestClassSetup)
  function fixture(test)
   [test.Eph,test.Config,test.Schedule]=v3SyntheticFixture('separated_visits');
   test.Features=ctocscreen.v3ScheduleFeatures(test.Schedule,test.Eph,test.Config);
  end
 end
 methods(Test)
  function evaporationPreservesPositiveExperience(test)
   c=test.Config; ph=newTable();
   [~,ph]=ctocscreen.v3Pheromone('deposit',ph,'success',.1,c);
   for k=1:100
    [~,ph]=ctocscreen.v3Pheromone('evaporate',ph,'',0,c);
   end
   unknown=ctocscreen.v3Pheromone('get',ph,'unknown',0,c);
   test.verifyGreaterThan(ph('success'),unknown);
   [~,ph]=ctocscreen.v3Pheromone('deposit',ph,'success',1e6,c);
   test.verifyEqual(ph('success'),c.pheromone_ceiling);
  end
  function repeatedAndNearbySchedulesDoNotMultiplyReward(test)
   f=test.Features; c=test.Config;
   [ph,ledger,first]=ctocscreen.v3Reinforce(newTable(),{},f,{'arc','arc'},10,c);
   [ph,ledger,again]=ctocscreen.v3Reinforce(ph,ledger,f,{'arc'},10,c);
   near=shiftFeature(f,.1,c);
   [ph,ledger,nearby]=ctocscreen.v3Reinforce(ph,ledger,near,{'arc'},11,c);
   test.verifyEqual(again.amount_per_key,0); test.verifyEqual(nearby.amount_per_key,0);
   test.verifyEqual(ph('arc'),c.pheromone_baseline+first.amount_per_key,'AbsTol',1e-14);
   test.verifyEqual(numel(ledger),1); test.verifyEqual(first.key_count,1);
   [ph,ledger,improved]=ctocscreen.v3Reinforce(ph,ledger,near,{'arc'},9,c);
   test.verifyEqual(improved.reason,'incremental_cost_gain');
   test.verifyEqual(ph('arc'),c.pheromone_baseline+c.pheromone_gain/(1+9),'AbsTol',1e-14);
   test.verifyEqual(ledger{1}.best_J,9);
  end
  function distinctFamilyEarnsFeedback(test)
   f=test.Features; c=test.Config;
   [ph,ledger]=ctocscreen.v3Reinforce(newTable(),{},f,{'arc'},10,c);
   [~,ledger,report]=ctocscreen.v3Reinforce(ph,ledger,shiftFeature(f,3,c),{'other'},12,c);
   test.verifyEqual(numel(ledger),2); test.verifyGreaterThan(report.amount_per_key,0);
  end
  function overlappingFamiliesUseBestKnownCost(test)
   f=test.Features; c=test.Config;
   [ph,ledger]=ctocscreen.v3Reinforce(newTable(),{},f,{'arc'},10,c);
   [ph,ledger]=ctocscreen.v3Reinforce(ph,ledger,shiftFeature(f,1.1,c),{'arc'},12,c);
   [~,~,report]=ctocscreen.v3Reinforce(ph,ledger,shiftFeature(f,.8,c),{'arc'},11,c);
   test.verifyEqual(report.amount_per_key,0);
  end
  function delayedBurnExperienceUsesActualDeparture(test)
   c=test.Config; c.root_count=1; c.actions_per_node=2; c.action_durations_s=300; c.budget_s=10;
   stream=RandStream('mt19937ar','Seed',15); roots=ctocscreen.v3Roots(test.Eph,c,stream);
   node=roots{1}; node.schedule=test.Schedule; node.schedule.witness_times_s(:)=NaN;
   node.state=test.Features.initial_state; node.visited(:)=false;
   children=ctocscreen.v3Expand(node,test.Eph,c,stream,newTable(),tic);
   burns=children(cellfun(@(n)~isempty(n.schedule.maneuver_times_s),children));
   test.assertNotEmpty(burns);
   child=burns{1}; s=child.schedule; tb=s.maneuver_times_s(1);
   test.verifyGreaterThan(tb,0);
   [x,~,arc]=ctocscreen.v3Arc(node.state,0,tb,test.Eph.model,c);
   scan=ctocscreen.v3ScanArc(arc,@(ids,t,mode)ctocscreen.v3QueryTargets(test.Eph,ids,t,mode),c,nan(35,1));
   actual=struct('state',x.','t',tb,'visited',scan.visited);
   expected=ctocscreen.v3StateKey(actual,2,0,s.duration_s-tb,c);
   test.verifyEqual(child.keys{end},expected);
  end
  function physicalFeaturesIgnoreAngleWrapAndWitnessRelabel(test)
   s=test.Schedule; s.initial_q(6)=s.initial_q(6)+2*pi;
   f=ctocscreen.v3ScheduleFeatures(s,test.Eph,test.Config);
   test.verifyLessThan(ctocscreen.v3FeatureDistance(f,test.Features,test.Config),1e-7);
   s=test.Schedule; s.witness_times_s=flipud(s.witness_times_s);
   f=ctocscreen.v3ScheduleFeatures(s,test.Eph,test.Config);
   test.verifyEqual(ctocscreen.v3FeatureDistance(f,test.Features,test.Config),0);
  end
  function topologyAndPathRemainDistinct(test)
   f=test.Features; b=f; b.physical_key='synthetic_path'; b.states(:,1)=b.states(:,1)+1000;
   test.verifyGreaterThanOrEqual(ctocscreen.v3FeatureDistance(f,b,test.Config),1);
   b=f; b.physical_key='synthetic_burn'; b.maneuver_times_s=200; b.delta_v_km_s=[.01 0 0];
   test.verifyGreaterThanOrEqual(ctocscreen.v3FeatureDistance(f,b,test.Config),1);
   b=shiftFeature(f,.1,test.Config); b.witness_times_s=flipud(b.witness_times_s);
   test.verifyGreaterThanOrEqual(ctocscreen.v3FeatureDistance(f,b,test.Config),1);
  end
  function nearCopyFloodCannotFillPool(test)
   c=test.Config; pool={work(test.Features,10)};
   for k=1:20
    [pool,report]=ctocscreen.v3SelectWorkPool(pool,work(shiftFeature(test.Features,k/100,c),11),c);
    test.verifyFalse(report.retained);
   end
   test.verifyEqual(numel(pool),1);
   [pool,report]=ctocscreen.v3SelectWorkPool(pool,work(shiftFeature(test.Features,.2,c),9),c);
   test.verifyTrue(report.retained); test.verifyEqual(numel(pool),1); test.verifyEqual(pool{1}.J,9);
  end
  function saUphillReplacesLocalWorkWithoutChangingBest(test)
   c=test.Config; best=work(test.Features,10); pool={best};
   [accepted,~,~]=ctocscreen.v3SAAccept(10,11,100,RandStream('mt19937ar','Seed',1));
   test.assertTrue(accepted);
   [pool,report]=ctocscreen.v3SelectWorkPool(pool,work(shiftFeature(test.Features,.1,c),11),c,true);
   test.verifyTrue(report.retained); test.verifyEqual(numel(pool),1);
   test.verifyEqual(pool{1}.J,11); test.verifyEqual(best.J,10);
  end
  function capacityKeepsBestAndSeparatedGeometry(test)
   c=test.Config; c.archive_size=2; f=test.Features;
   pool={work(f,1),work(shiftFeature(f,2,c),2)};
   [pool,report]=ctocscreen.v3SelectWorkPool(pool,work(shiftFeature(f,10,c),3),c);
   test.verifyTrue(report.retained); test.verifyEqual(cellfun(@(w)w.J,pool),[1 3]);
   [pool,report]=ctocscreen.v3SelectWorkPool(pool,work(shiftFeature(f,5,c),4),c,true);
   test.verifyTrue(report.retained); test.verifyEqual(sort(cellfun(@(w)w.J,pool)),[1 4]);
   c.archive_size=1;
   [pool,report]=ctocscreen.v3SelectWorkPool(pool,work(shiftFeature(f,20,c),6),c,true);
   test.verifyTrue(report.retained); test.verifyEqual(numel(pool),1); test.verifyEqual(pool{1}.J,6);
  end
  function stateKeysSeparateDirectionsAndVisitSets(test)
   a=struct('state',[7000 0 0 0 7 0],'t',0,'visited',false(35,1));
   b=a; b.state=[0 7000 0 -7 0 0];
   key=@(n)ctocscreen.v3StateKey(n,3,35,600,test.Config);
   test.verifyNotEqual(key(a),key(b));
   a.visited([1 4])=true; b=a; b.visited(:)=false; b.visited([2 3])=true;
   test.verifyNotEqual(key(a),key(b));
  end
  function experienceCreditsEveryVisitAndCoast(test)
   s=test.Schedule; s.maneuver_times_s=25; s.delta_v_km_s=[1e-5 0 0];
   checked=ctocscreen.v3Verify(s,test.Eph,test.Config); test.assertTrue(checked.passed);
   s.witness_times_s=checked.witness_times_s;
   [keys,features]=ctocscreen.v3ScheduleExperience(s,test.Eph,test.Config);
   [~,trace]=ctocscreen.v3Replay(s,test.Eph,test.Config,false);
   node=struct('state',trace.states(2,:),'t',25,'visited',s.witness_times_s<=25);
   node.state(4:6)=node.state(4:6)-s.delta_v_km_s;
   for target=1:35
    expected=ctocscreen.v3StateKey(node,3,target,s.witness_times_s(target)-25,test.Config);
    test.verifyTrue(ismember(expected,keys));
   end
   test.verifyTrue(ismember(ctocscreen.v3StateKey(node,2,0,575,test.Config),keys));
   node.state=trace.states(1,:); node.t=0; node.visited(:)=false;
   test.verifyTrue(ismember(ctocscreen.v3StateKey(node,1,0,25,test.Config),keys));
   test.verifyEqual(features.physical_key,ctocscreen.v3ScheduleKey(s));
  end
  function beamMergesCopiesBeforeQuotas(test)
   a=struct('state',[7000 0 0 0 7 0],'t',0,'visited',false(35,1),'estimate',1,'last_gain',1);
   a.visited(1:5)=true;
   coast=a; coast.state=[0 7000 0 -7 0 0]; coast.last_gain=0; coast.estimate=100;
   diverse=a; diverse.state=-a.state; diverse.estimate=200;
   chosen=ctocscreen.v3SelectBeam([repmat({a},1,20),{coast,diverse}],4,test.Config);
   test.verifyEqual(numel(chosen),3);
   test.verifyTrue(any(cellfun(@(n)n.last_gain==0,chosen)));
   test.verifyTrue(any(cellfun(@(n)isequal(n.state,diverse.state),chosen)));
   chosen=ctocscreen.v3SelectBeam({coast,a},1,test.Config);
   test.verifyEqual(numel(chosen),1); test.verifyEqual(chosen{1}.estimate,1);
  end
  function searchGatesAndCoalescesInitialSeeds(test)
   c=test.Config; c.budget_s=10; c.root_count=1; c.colony_count=1;
   near=test.Schedule; near.initial_q(6)=near.initial_q(6)+1e-6;
   bad=test.Schedule; bad.initial_q(6)=bad.initial_q(6)+.01;
   c.initial_candidates={test.Schedule,near,bad};
   c.stop_file=mfilename('fullpath'); c.stop_file=[c.stop_file '.m'];
   sim=fileparts(fileparts(mfilename('fullpath')));
   folder=tempname(fullfile(sim,'runs','v3','development'));
   state=ctocscreen.v3Search(folder,test.Eph,c);
   test.verifyEqual(state.iteration,0); test.verifyEqual(numel(state.initial_seed_checks),3);
   test.verifyTrue(state.initial_seed_checks{1}.passed); test.verifyTrue(state.initial_seed_checks{2}.passed);
   test.verifyFalse(state.initial_seed_checks{3}.passed);
   test.verifyEqual(numel(state.work_pool),1); test.verifyEqual(numel(state.reinforcement_ledger),1);
   test.verifyEqual(numel(state.feedback_history),2);
   test.verifyEqual(state.feedback_history{2}.amount_per_key,0);
   test.verifyEqual(state.elite.verification.total_dv_km_s,0);
   test.verifyEqual(state.signature.schema,'free_maneuver_search_v3_8');
  end
 end
end

function ph=newTable()
ph=containers.Map('KeyType','char','ValueType','double');
end

function b=shiftFeature(a,scale,c)
b=a; b.physical_key=sprintf('synthetic_shift_%g',scale);
b.initial_state(1)=b.initial_state(1)+scale*c.diversity_scales(1);
end

function w=work(f,J)
w=struct('schedule',struct('synthetic_test_only',true),'J',J,'features',f);
end
