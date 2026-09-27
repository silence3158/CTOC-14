function report=check_v4_incumbent_tools(label)
%CHECK_V4_INCUMBENT_TOOLS Targeted checks for the budget filter and re-expansion.
% Cold nodes only (seed 888 roots and short greedy chains built here); no
% historical trajectory is read. Questions:
%  (a) With a tight remaining budget R, does max_leg_dv=R+margin turn wasted
%      over-budget guidance into admissible children (same node, same RNG)?
%  (b) Does a second expansion with exclude=tried give new encounters only?
%  (c) Scheduler: with return limit 4, landings 1-3 keep the node, 4 stops it.
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'),fullfile(sim,'scripts'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(); stream=RandStream('mt19937ar','Seed',888); margin=0.1;
% (c) synthetic scheduler scenario, limit 4 (burn counts 0..8 along a chain).
R=zeros(1,12); S=false(1,12);
for landing=1:3
 [p,R,S,~,h]=v4_bound4_return_landings(1:9,0:8,R,S,4,4);
 assert(isequal(p,1:5)&&R(5)==landing&&~S(5)&&isempty(h),'Landing %d must keep the node.',landing);
end
[p,R,S,legs,h]=v4_bound4_return_landings(1:9,0:8,R,S,4,4);
assert(isequal(p,1)&&R(5)==4&&S(5)&&isequal(h,5)&&size(legs,1)==2&&R(1)==1,'Fourth landing must stop and cascade.');
fprintf('SCHEDULER limit 4: landings 1-3 keep, 4th stops and cascades\n');
% Cold nodes: four roots, two greedy chains of six cheapest-child steps.
nodes={}; memory=[];
for k=1:4, nodes{end+1}=ctocscreen.v4.root(k,eph,c,stream); end %#ok<AGROW>
for k=1:2
 n=nodes{k};
 for step=1:6
  [kids,~,memory]=ctocscreen.v4.expand(n,eph,c,stream,memory,c.action_seconds);
  kids=kids(cellfun(@(x)x.actual.visit_count>n.actual.visit_count,kids));
  if isempty(kids), break; end
  [~,j]=min(cellfun(@(x)x.actual.total_dv_km_s,kids)); n=kids{j};
  if step==3||step==6, nodes{end+1}=n; end %#ok<AGROW>
 end
end
budgets=[0.3 0.6 1.2]; rows=zeros(0,9); ex=zeros(0,5);
for i=1:numel(nodes)
 n=nodes{i}; J0=n.actual.total_dv_km_s;
 for b=budgets
  s0=stream.State;
  t1=tic; [kidsA,rA]=ctocscreen.v4.expand(n,eph,c,stream,[],c.action_seconds); tA=toc(t1);
  stream.State=s0;
  t2=tic; [kidsB,rB]=ctocscreen.v4.expand(n,eph,c,stream,[],c.action_seconds,struct('max_leg_dv',b+margin)); tB=toc(t2);
  % Only children that add a visit count; a zero-cost coast fallback does not.
  [legA,gainA]=legCosts(kidsA,J0,n.actual.visit_count); [legB,gainB]=legCosts(kidsB,J0,n.actual.visit_count);
  rows(end+1,:)=[i,n.actual.visit_count,b,sum(gainA&legA<=b),sum(gainA&legA>b), ...
   sum(gainB&legB<=b),sum(gainB&legB>b),rB.budget_filtered,tB-tA]; %#ok<AGROW>
 end
 % (b) re-expansion with the first call's tried encounters excluded.
 [k1,r1]=ctocscreen.v4.expand(n,eph,c,stream,[],c.action_seconds);
 [k2,r2]=ctocscreen.v4.expand(n,eph,c,stream,[],c.action_seconds,struct('exclude',r1.tried));
 clash=0;
 for j=1:size(r2.tried,1)
  clash=clash+any(r1.tried(:,1)==r2.tried(j,1)&abs(r1.tried(:,2)-r2.tried(j,2))<600);
 end
 assert(clash==0,'Re-expansion retried an excluded encounter.');
 k1=k1(cellfun(@(x)x.actual.visit_count>n.actual.visit_count,k1));
 k2=k2(cellfun(@(x)x.actual.visit_count>n.actual.visit_count,k2));
 newKeys=setdiff(cellfun(@(x)ctocscreen.v4.controlKey(x.q),k2,'UniformOutput',false), ...
  cellfun(@(x)ctocscreen.v4.controlKey(x.q),k1,'UniformOutput',false));
 ex(end+1,:)=[i,numel(k1),numel(k2),numel(newKeys),r2.excluded]; %#ok<AGROW>
end
fprintf('FILTER per call (node visits, R): admissible without/with filter, over-budget without/with, filtered proposals\n');
for k=1:size(rows,1)
 fprintf('  node %d (%2d visits) R=%.1f: admissible %d -> %d, over-budget %d -> %d, filtered %d\n',rows(k,[1 2 3 4 6 5 7 8]));
end
fprintf('FILTER totals: admissible %d -> %d, over-budget %d -> %d over %d calls\n', ...
 sum(rows(:,4)),sum(rows(:,6)),sum(rows(:,5)),sum(rows(:,7)),size(rows,1));
fprintf('REEXPAND visit-adding children first/second, new controls, excluded proposals:\n'); disp(ex);
report=struct('diagnostic_only',true,'cold_nodes_only',true,'filter_rows',rows, ...
 'filter_columns',{{'node','visits','R','adm_no_filter','over_no_filter','adm_filter','over_filter','filtered','dt_s'}}, ...
 're_expand',ex,'margin_km_s',margin,'signature',ctocscreen.v4.signature());
save(fullfile(sim,'runs/v4/development',['incumbent_tools_' label '.mat']),'report');
fprintf('V4 INCUMBENT TOOL CHECKS PASSED\n');
end
function [leg,gain]=legCosts(kids,J0,visits)
leg=cellfun(@(x)x.actual.total_dv_km_s-J0,kids);
gain=cellfun(@(x)x.actual.visit_count>visits,kids);
end
