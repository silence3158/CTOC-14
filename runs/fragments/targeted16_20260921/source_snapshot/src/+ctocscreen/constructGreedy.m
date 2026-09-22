function [c,diagnostic]=constructGreedy(p,cfg,stream)
%CONSTRUCTGREEDY Randomized restricted-candidate construction; not an optimum.
base=ctocscreen.constructCandidates(p,struct(),stream,1);
seedTarget=[];
if isfield(cfg,'geometry_seed_probability') && rand(stream)<cfg.geometry_seed_probability
 [base,seedTarget]=ctocscreen.geometrySeed(base,p,stream);
end
s=ctocscreen.initialState(base.initial_q,p.mu_km3_s2,p.re_km);
[s,~]=ctocscreen.propagateTwoBody(s,base.wait_s,p.mu_km3_s2);
t=base.wait_s; remaining=1:35; order=zeros(1,35); dt=zeros(1,35);
refs=cell(1,35); diagnostic=struct('status','ok','failure_leg',0); c=[];
for k=1:35
 maxTime=p.horizon_s-t-(35-k)*cfg.min_tof_s;
 if isfield(cfg,'construct_time_budget_factor')
  maxTime=min(maxTime,cfg.construct_time_budget_factor*(p.horizon_s-t)/(36-k));
 end
 idx=remaining(randperm(stream,numel(remaining)));
 trial=idx(1:min(cfg.construct_target_count,numel(idx)));
 if k==1 && ~isempty(seedTarget), trial=seedTarget; end
 options=[]; records={};
 % If the restricted target set fails, retry all remaining targets.
 for pass=1:2
  for id=trial
   times=cfg.construct_times_s.*exp(0.2*randn(stream,1,numel(cfg.construct_times_s)));
   times=min(times,maxTime);
   times=unique(times(times>=cfg.min_tof_s));
   for tof=times
    [target,~]=ctocscreen.targetStates(p,id,t+tof);
    [bs,~]=ctocscreen.enumerateBranches(s(1:3),target(1:3),tof,p.mu_km3_s2,cfg.branch_policy);
    for j=1:numel(bs)
     a=ctocscreen.checkArc([s(1:3) bs(j).v_depart],tof,p);
     floorAltitude=200;
     if k>1, floorAltitude=max(200,cfg.construction_min_altitude_after_first_km); end
     if strcmp(a.status,'ok') && a.min_altitude_km>=floorAltitude
      options(end+1)=norm(bs(j).v_depart-s(4:6)); %#ok<AGROW>
      records{end+1}={id,tof,bs(j)}; %#ok<AGROW>
     end
    end
   end
  end
  if ~isempty(options), break; end
  trial=remaining;
 end
 if isempty(options), diagnostic.status='no_branch'; diagnostic.failure_leg=k; return; end
 if isfield(cfg,'construct_refine_count') && cfg.construct_refine_count>0
  [~,rank]=sort(options); refinements=rank(1:min(cfg.construct_refine_count,numel(rank)));
  for index=refinements
   rec=records{index}; id=rec{1}; anchor=rec{2}; ref=rec{3};
   lower=max(cfg.min_tof_s,anchor*.65);
   upper=min(anchor*1.4,maxTime);
   if upper<=lower, continue; end
   [optTime,optCost]=fminbnd(@timeCost,lower,upper,optimset('Display','off','TolX',.5,'MaxFunEvals',24));
   [checkedCost,optBranch]=timeCost(optTime);
   if ~isempty(optBranch) && optCost<options(index) && checkedCost<options(index)
    options(index)=checkedCost; records{index}={id,optTime,optBranch};
   end
  end
 end
 [~,ii]=sort(options); pick=ii(1);
 if rand(stream)<0.2, pick=ii(randi(stream,min(3,numel(ii)))); end
 rec=records{pick}; order(k)=rec{1}; dt(k)=rec{2}; refs{k}=rec{3};
 [s,info]=ctocscreen.propagateTwoBody([s(1:3) refs{k}.v_depart],dt(k),p.mu_km3_s2);
 if ~strcmp(info.status,'ok'), diagnostic.status='solver_failure'; diagnostic.failure_leg=k; return; end
 t=t+dt(k); remaining(remaining==order(k))=[];
end
c=ctocscreen.makeCandidate(base.initial_q,order,base.wait_s,dt, ...
 cellfun(@(x)x.branch_id,refs),base.seed,'','randomized_greedy');
c.branch_refs=refs;
 function [cost,branch]=timeCost(tof)
  cost=1e6; branch=[];
  target=ctocscreen.targetStates(p,id,t+tof);
  bs=ctocscreen.enumerateBranches(s(1:3),target(1:3),tof,p.mu_km3_s2,cfg.branch_policy);
  jj=find([bs.branch_id]==ref.branch_id,1);
  if isempty(jj), return; end
  arc=ctocscreen.checkArc([s(1:3) bs(jj).v_depart],tof,p);
  floorAlt=200; if k>1, floorAlt=max(200,cfg.construction_min_altitude_after_first_km); end
  if ~strcmp(arc.status,'ok') || arc.min_altitude_km<floorAlt, return; end
  branch=bs(jj); cost=norm(branch.v_depart-s(4:6));
 end
end
