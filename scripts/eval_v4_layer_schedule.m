function report=eval_v4_layer_schedule(seed,budget,width,label)
%EVAL_V4_LAYER_SCHEDULE Evaluation of a proposed schedule, not a search entry.
% Layered-beam construction only: every round expands every beam member
% (expand.m), then selectBeam keeps `width` nodes. No B, no shared arcs,
% no restructuring. Cold start, empty history. Every complete candidate is
% independently replayed; the report states actual independent results.
if nargin<2, budget=300; end
if nargin<3, width=16; end
if nargin<4, label=sprintf('s%d',seed); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(struct('seed',seed,'beam_width',width,'budget_s',budget));
stream=RandStream('mt19937ar','Seed',seed); clock=tic; memory=[];
beam={}; nextRoot=1;
for k=1:c.root_count, beam{end+1}=ctocscreen.v4.root(nextRoot,eph,c,stream); nextRoot=nextRoot+1; end %#ok<AGROW>
deadline=budget-c.verify_reserve_s; rounds=0; complete={}; history=zeros(0,4); verified=[];
while toc(clock)<deadline&&~isempty(beam)
 rounds=rounds+1; next={};
 if mod(rounds,c.root_every)==0
  beam{end+1}=ctocscreen.v4.root(nextRoot,eph,c,stream); nextRoot=nextRoot+1; %#ok<AGROW>
 end
 for k=1:numel(beam)
  if toc(clock)>=deadline, break; end
  n=beam{k}; if n.actual.visit_count>=35||n.q.T>=eph.model.horizon_s-1, continue; end
  [children,~,memory]=ctocscreen.v4.expand(n,eph,c,stream,memory,min(c.action_seconds,deadline-toc(clock)));
  for j=1:numel(children)
   if children{j}.actual.passed, complete{end+1}=children{j}; end %#ok<AGROW>
  end
  next=[next,children]; %#ok<AGROW>
 end
 if isempty(next), break; end
 [beam,~]=ctocscreen.v4.selectBeam(next,c,eph);
 counts=cellfun(@(n)n.actual.visit_count,beam); [kmax,j]=max(counts);
 history(end+1,:)=[toc(clock),kmax,beam{j}.actual.total_dv_km_s,numel(complete)]; %#ok<AGROW>
end
searchSeconds=toc(clock); firstComplete=NaN;
if ~isempty(complete)
 costs=cellfun(@(n)n.actual.total_dv_km_s,complete); [~,j]=min(costs); best=complete{j};
 firstComplete=history(find(history(:,4)>0,1),1);
else
 counts=cellfun(@(n)n.actual.visit_count,beam); costs=cellfun(@(n)n.actual.total_dv_km_s,beam);
 [~,order]=sortrows([-counts(:),costs(:)]); best=beam{order(1)};
end
[verified,~]=ctocscreen.v4.replay(best.q,eph,c,true);
fprintf('LAYER seed=%d width=%d rounds=%d search=%.1f s complete_found=%d first_complete=%.1f s\n', ...
 seed,width,rounds,searchSeconds,numel(complete),firstComplete);
fprintf('LAYER independent %d/35 dv=%.9f passed=%d height=%.3f total=%.1f s\n', ...
 verified.visit_count,verified.total_dv_km_s,verified.passed,verified.min_altitude_lower_km,toc(clock));
thin=@(n)rmfield(n,'trace');
report=struct('evaluation_only',true,'seed',seed,'width',width,'budget_s',budget,'rounds',rounds, ...
 'search_seconds',searchSeconds,'history',history,'complete_found',numel(complete), ...
 'first_complete_s',firstComplete,'best',thin(best),'verification',verified, ...
 'complete_q',{cellfun(@(n)n.q,complete,'UniformOutput',false)},'signature',ctocscreen.v4.signature());
out=fullfile(sim,'runs/v4/development/schedule_eval'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,['layer_' label '.mat']),'report');
end
