function report=check_v4_core()
%CHECK_V4_CORE Targeted numerical contracts; synthetic checks are not scores.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'),fullfile(sim,'tests'));
[eph,~,old]=v3SyntheticFixture('separated_visits');
c=ctocscreen.v4.defaults(struct('max_step_s',30,'scan_step_s',60,'joint_iterations',8,'solver_seconds',3));
q=struct('x0',ctocscreen.initialState(old.initial_q,eph.model.mu,eph.model.re).', ...
 'tau',[90;360],'u',[.002 -.001 .0005;-.0004 .0003 .0002],'T',600,'witness',nan(35,1));
[a,tr]=ctocscreen.v4.replay(q,eph,c); assert(~strcmp(a.status,'propagation_failure'),a.failure_reason);
node=struct('q',q,'actual',a,'trace',tr,'root_id',1,'root_kind','synthetic', ...
 'seed_target',1,'seed_duration',100,'attempts',0,'zero_gain',0,'generation',0,'origin','synthetic','heuristic_H',NaN);
ids=[8;24]; theta=[160;420]; p=ctocscreen.v4.problem(node,ids,theta,'full',300,eph,c);
z=p.z0; e=ctocscreen.v4.evaluate(p,z,true);
assert(e.max_position_defect<1e-5,'Initial nodes must follow fixed controls.');
rng(888); direction=randn(p.n,1); direction=direction/norm(direction); h=2e-7;
% Keep the initial epoch/end bounds and exact fixed dimensions untouched.
direction(p.it(end))=0; direction(~p.free)=0;
plus=ctocscreen.v4.evaluate(p,z+h*direction,false); minus=ctocscreen.v4.evaluate(p,z-h*direction,false);
numeric=[(plus.eq-minus.eq);(plus.g-minus.g);(plus.visit(:)-minus.visit(:));(plus.ecc-minus.ecc)]/(2*h);
analytic=[e.Jeq;e.Jg;e.Jvisit;e.Jecc]*direction;
relative=norm(numeric-analytic,Inf)/max(1,norm(numeric,Inf));
assert(relative<3e-5,'Analytic directional Jacobian mismatch: %.3g',relative);
fprintf('V4 CHECK Jacobian relative error %.3g\n',relative);
local=ctocscreen.v4.problem(node,ids,theta,'local',450,eph,c);
assert(~any(local.free(local.iu(:,1))),'Earlier pulse should be frozen in this local task.');
assert(all(local.free(local.ix(:,end))),'Suffix state must not be frozen with controls.');
changed=q; changed.u(1,1)=changed.u(1,1)+1e-4;
[a2,tr2]=ctocscreen.v4.replay(changed,eph,c);
assert(norm(tr2.pre(:,end)-tr.pre(:,end))>.001,'Earlier controls must change the suffix.');
assert(~strcmp(ctocscreen.v4.controlKey(changed),ctocscreen.v4.controlKey(q)));
zero=q; zero.tau=[zero.tau;240]; zero.u=[zero.u;zeros(1,3)];
[az,tz]=ctocscreen.v4.replay(zero,eph,c);
assert(norm(tz.post(:,end)-tr.post(:,end))<1e-8&&abs(az.total_dv_km_s-a.total_dv_km_s)<1e-12);
coincident=ctocscreen.v4.problem(node,8,90,'full',90,eph,c);
ec=ctocscreen.v4.evaluate(coincident,coincident.z0,true);
assert(ec.max_position_defect<1e-5,'Coincident visit/burn event defect.');
fprintf('V4 CHECK event identity and changed suffix passed\n');
step=ctocscreen.v4.coneStep(p,z,e,.03,10,false,3);
assert(step.ok,'SOCP did not produce a feasible linearized step.');
[pred,~,~]=ctocscreen.v4.modelMerit(p,z,e,step.d,10,false);
assert(pred<e.V,'Restoration subproblem did not lower predicted violation.');
fprintf('V4 CHECK cone exit=%d time=%.3f predicted %.6g -> %.6g\n',step.exitflag,step.seconds,e.V,pred);
[out,jr]=ctocscreen.v4.joint(node,ids,theta,'full',300,eph,c,15);
assert(jr.iterations>0&&~strcmp(jr.status,'numerical_failure'),jr.failure_reason);
assert(~isempty(out)&&out.actual.independent==false);
assert(strcmp(out.actual.dataset_kind,'synthetic_test_only'));
fprintf('V4 CHECK joint iterations=%d accepted=%d active_pass=%d status=%s\n',jr.iterations,jr.accepted_steps,jr.active_passed,jr.status);
% Prefix reuse is compared to a fresh replay of identical controls.
short=q; short.T=450; short.tau=q.tau; short.u=q.u;
[sa,st]=ctocscreen.v4.replay(short,eph,c); parent=node; parent.q=short; parent.actual=sa; parent.trace=st;
long=short; long.T=600; long.tau(end+1,1)=450; long.u(end+1,:)=[.0002 0 0];
[cached,ct]=ctocscreen.v4.replay(long,eph,c,false,parent); [fresh,ft]=ctocscreen.v4.replay(long,eph,c);
reuseError=norm(ct.post(:,end)-ft.post(:,end));
assert(cached.reused_prefix&&reuseError<1e-5&&cached.visit_count==fresh.visit_count);
long.tau(end)=480; [delayed,dt]=ctocscreen.v4.replay(long,eph,c,false,parent);
[fromZero,zt]=ctocscreen.v4.replay(long,eph,c);
% Since 2026-09-26 reuse may stop at the parent's last impulse; it must still
% equal the epoch replay exactly.
assert(isequal(dt.post,zt.post)&&delayed.visit_count==fromZero.visit_count);
% A preceding screen must not suppress independent positive reinforcement.
reward=node; reward.q.u=q.u/100; [reward.actual,reward.trace]=ctocscreen.v4.replay(reward.q,eph,c);
reward.actual.passed=true; reward.actual.independent=false;
[mem,~]=ctocscreen.v4.feedback([],'observe',reward,c);
reward.actual.independent=true; [mem,~]=ctocscreen.v4.feedback(mem,'observe',reward,c);
assert(mem.positive_updates==1&&any(mem.positive>0));
[mem,~]=ctocscreen.v4.feedback(mem,'observe',reward,c); assert(mem.positive_updates==1);
report=struct('dataset_kind','synthetic_test_only','passed',true,'jacobian_relative_error',relative, ...
 'prefix_reuse_error',reuseError,'cone_exitflag',step.exitflag,'cone_seconds',step.seconds,'joint',jr);
folder=fullfile(sim,'runs','v4','development'); if ~isfolder(folder), mkdir(folder); end
save(fullfile(folder,'core_checks.mat'),'report');
fprintf('V4 CORE CHECKS PASSED (synthetic only)\n');
end
