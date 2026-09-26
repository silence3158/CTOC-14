function report=check_v4_tail(label)
%CHECK_V4_TAIL Targeted checks for the tail-window B task.
% Q1 Jacobian of the tail transcription matches finite differences.
% Q2 decode(z0) reproduces the original controls exactly.
% Q3 On real complete/partial chains (diagnostic inputs only), does tail B
%    produce replay-accepted Delta-V savings, where full-history B did not (E3)?
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults();
files={'runs/v4/development/crossing_t2_w16_w16_defaults.mat','runs/v4/development/crossing_t2_w8_w8_defaults.mat'};
rows=zeros(0,8);
for f=1:numel(files)
 s=load(fullfile(sim,files{f})); q=s.report.q;
 [a,tr]=ctocscreen.v4.replay(q,eph,c);
 node=struct('q',tr.q,'actual',a,'trace',tr,'root_id',1,'root_kind','e','seed_target',1,'seed_duration',3600, ...
  'attempts',0,'zero_gain',0,'generation',0,'origin','e','heuristic_H',NaN);
 node.q.witness=a.witness_times_s; ids=find(a.distance_km<=1); M=numel(node.q.tau);
 if f==1
  p=ctocscreen.v4.problem(node,ids,a.witness_times_s(ids),'tail',node.q.tau(M-2),eph,c);
  q0=ctocscreen.v4.decode(p,p.z0);
  assert(isequal(q0.u,node.q.u)&&isequal(q0.tau,node.q.tau)&&isequal(q0.x0,node.q.x0)&&q0.T==node.q.T,'Tail decode changed controls.');
  e=ctocscreen.v4.evaluate(p,p.z0,true);
  fprintf('TAIL chart: N=%d M=%d K=%d n=%d t0=%.1f gap_km=%.3g\n',p.N,p.M,p.K,p.n,p.t0,e.model_gap_km);
  % Arrival witnesses coincide with the next burn or T; a two-sided difference
  % crosses the ordering clamp there, so only noncoincident times are perturbed.
  tt=p.raw_times; co=false(size(tt));
  for k=1:numel(tt), co(k)=any(abs(tt([1:k-1,k+1:end])-tt(k))<1e-6); end
  rng(3); d=randn(p.n,1); d(~p.free)=0; d(p.it(co))=0; d(p.it(end))=0; d=d/norm(d); h=2e-7;
  ep=ctocscreen.v4.evaluate(p,p.z0+h*d,false); em=ctocscreen.v4.evaluate(p,p.z0-h*d,false);
  num=[(ep.eq-em.eq);(ep.g-em.g);(ep.visit(:)-em.visit(:))]/(2*h); ana=[e.Jeq;e.Jg;e.Jvisit]*d;
  rel=norm(num-ana,Inf)/max(1,norm(num,Inf)); assert(rel<3e-5,'Tail Jacobian mismatch %.3g',rel);
  fprintf('TAIL Jacobian relative error %.3g\n',rel);
 end
 for w=[2 4]
  for budget=[8 30]
   t1=tic; [child,jr]=ctocscreen.v4.joint(node,ids,a.witness_times_s(ids),'tail',node.q.tau(M-w+1),eph,c,budget);
   J1=a.total_dv_km_s; if ~isempty(child)&&jr.active_passed, J1=child.actual.total_dv_km_s; end
   rows(end+1,:)=[f,w,budget,a.total_dv_km_s,J1,jr.iterations,jr.accepted_steps,toc(t1)]; %#ok<AGROW>
   fprintf('TAIL %s last %d burns budget %d: %.6f -> %.6f (status %s, it %d, acc %d, %.1f s)\n', ...
    files{f},w,budget,a.total_dv_km_s,J1,jr.status,jr.iterations,jr.accepted_steps,toc(t1));
  end
 end
end
% Independent confirmation of the best saving found.
[~,j]=max(rows(:,4)-rows(:,5));
report=struct('diagnostic_only',true,'rows',rows,'signature',ctocscreen.v4.signature());
if rows(j,4)-rows(j,5)>1e-7
 s=load(fullfile(sim,files{rows(j,1)})); q=s.report.q; [a,tr]=ctocscreen.v4.replay(q,eph,c);
 node=struct('q',tr.q,'actual',a,'trace',tr,'root_id',1,'root_kind','e','seed_target',1,'seed_duration',3600, ...
  'attempts',0,'zero_gain',0,'generation',0,'origin','e','heuristic_H',NaN); node.q.witness=a.witness_times_s;
 ids=find(a.distance_km<=1); M=numel(node.q.tau);
 child=ctocscreen.v4.joint(node,ids,a.witness_times_s(ids),'tail',node.q.tau(M-rows(j,2)+1),eph,c,rows(j,3));
 v=ctocscreen.v4.replay(child.q,eph,c,true);
 fprintf('TAIL independent check of best: %d/35 passed=%d dv=%.9f\n',v.visit_count,v.passed,v.total_dv_km_s);
 report.independent=v;
end
save(fullfile(sim,'runs/v4/development',['tail_checks_' label '.mat']),'report');
end
