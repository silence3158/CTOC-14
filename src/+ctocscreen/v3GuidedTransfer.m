function [dv,stats,alternatives]=v3GuidedTransfer(x,t0,t1,goal,m,c,clock,departureSeed)
%V3GUIDEDTRANSFER Retry branches, screen deep-Earth seeds, damp J2 correction.
x=x(:); goal=goal(:); alternatives={};
assert(toc(clock)<c.budget_s,'ctocscreen:v3:guidedBudget','Guided budget exhausted before enumeration.');
policy=struct('max_revolutions',ctocscreen.v3LambertRevolutions(x(1:3),goal,t1-t0,m.mu,c.guided_max_revolutions),'endpoint_tol_km',.001);
if nargin>=8&&~isempty(departureSeed)
 branches=struct('v_depart',departureSeed(:).');
else
 branches=ctocscreen.v3LambertBranches(x(1:3),goal,t1-t0,m.mu,policy);
end
assert(~isempty(branches),'ctocscreen:v3:noBranch','No Lambert seed.');
scores=arrayfun(@(b)norm(b.v_depart-x(4:6).'),branches);
for j=1:numel(branches)
 change=ctocscreen.v3PlaneChange(x,branches(j).v_depart-x(4:6).',c);
 scores(j)=scores(j)+c.plane_penalty_km_s*change.penalty;
end
[~,order]=sort(scores);
stats=struct('branches_tried',0,'screened_deep_earth',0,'last_error',''); dv=[];
for b=order(1:min(c.guided_branches,numel(order)))
 if toc(clock)>=c.budget_s, break; end
 stats.branches_tried=stats.branches_tried+1;
 try
  v=branches(b).v_depart(:);
  % Proposal screen only, NOT a J2 altitude certificate. Leave a wide margin.
  quick=ctocscreen.checkArc([x(1:3);v].',t1-t0,struct('mu_km3_s2',m.mu,'re_km',m.re));
  if strcmp(quick.status,'ok')&&quick.min_altitude_km<0
   stats.screened_deep_earth=stats.screened_deep_earth+1; continue
  end
  for iter=1:10
   if toc(clock)>=c.budget_s, break; end
   [y,P,sol]=ctocscreen.v3Arc([x(1:3);v],t0,t1,m,c,true);
   residual=y(1:3)-goal; distance=norm(residual);
   if distance<.05
    % Commit an accurate six-state arc, not the derivative integrator's state.
    % Full verification still replays every fixed pulse afresh from the epoch.
    [y,~,sol]=ctocscreen.v3Arc([x(1:3);v],t0,t1,m,c,false,true);
    residual=y(1:3)-goal; distance=norm(residual);
   end
   if distance<.05
    h=ctocscreen.v3Height(sol,m,c);
    if h.passed
     impulse=v-x(4:6);
     alternatives{end+1}=struct('dv',impulse,'final_state',y.','branch',b,'distance_km',distance, ...
      'arc',struct('solution',sol,'height',h));
     if isempty(dv), dv=impulse; end
     if nargout<3, return; end
    end
    break
   end
   if toc(clock)>=c.budget_s, break; end
   B=P(1:3,4:6); assert(rcond(B)>1e-12,'ctocscreen:v3:shootSingular','Singular correction.');
   step=B\residual; step=step*min(1,2/max(norm(step),eps)); accepted=false;
   for alpha=[1 .5 .25 .125]
    if toc(clock)>=c.budget_s, break; end
    try
     trial=ctocscreen.v3Arc([x(1:3);v-alpha*step],t0,t1,m,c);
     if norm(trial(1:3)-goal)<distance, v=v-alpha*step; accepted=true; break; end
    catch
     % Reject this trial step; the successful arc is still checked in J2.
    end
   end
   if ~accepted, break; end
  end
 catch err
  stats.last_error=err.identifier;
 end
end
stats.successful_branches=numel(alternatives);
if isempty(dv)&&nargin>=8&&~isempty(departureSeed)&&toc(clock)<c.budget_s
 % A time-ranked seed may fail the J2/height gate; retain other Lambert paths.
 [dv,retry,alternatives]=ctocscreen.v3GuidedTransfer(x,t0,t1,goal,m,c,clock);
 stats.branches_tried=stats.branches_tried+retry.branches_tried;
 stats.screened_deep_earth=stats.screened_deep_earth+retry.screened_deep_earth;
 stats.last_error=retry.last_error; stats.successful_branches=numel(alternatives);
end
if isempty(dv)
 error('ctocscreen:v3:guidedExhausted','Guided branches exhausted (%d tried, %d deep-Earth). %s', ...
 stats.branches_tried,stats.screened_deep_earth,stats.last_error);
end
end
