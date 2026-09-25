function report_v3_round_diag(label)
%REPORT_V3_ROUND_DIAG Read-only pipeline diagnostics for one V3 search run.
% Diagnostic only: it never feeds a search and never claims a solution.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
folder=fullfile(sim,'runs','v3','search',label);
S=load(fullfile(folder,'checkpoint.mat'),'state'); S=S.state;
fprintf('=== %s ===\n',label);
fprintf('seed %d | wall %.3f s | overrun %.3f | iterations %d\n',S.config.seed,S.elapsed_s,S.budget_overrun_s,S.iteration);
fprintf('sources: %d\n',numel(S.source_history));
outcomes=cellfun(@(h)h.outcome,S.source_history,'UniformOutput',false);
[u,~,idx]=unique(outcomes); cnt=accumarray(idx(:),1); [~,o]=sort(cnt,'descend');
for k=reshape(o,1,[])
 fprintf('  %s = %d\n',u{k},cnt(k));
end
% construction
if ~isempty(S.construction_history)
 ch=S.construction_history;
 fprintf('construction calls %d | total %.3f s | checks %.3f s over %d checks (%d cache hits)\n', ...
  numel(ch),sum(cellfun(@(x)x.elapsed_s,ch)),sum(cellfun(@(x)x.prefix_check_seconds,ch)), ...
  sum(cellfun(@(x)x.consistency_checks,ch)),sum(cellfun(@(x)x.consistency_cache_hits,ch)));
 fprintf('construction end visits: '); fprintf('%d ',cellfun(@(x)x.end_visits,ch)); fprintf('\n');
 fprintf('completion output counts: '); fprintf('%d ',cellfun(@(x)x.output_count,ch)); fprintf('\n');
end
% full construction (cold completion of a 35-target seed)
if ~isempty(S.full_construction_history)
 fh=S.full_construction_history;
 fprintf('full-construction calls %d | total %.3f s | end visits ',numel(fh),sum(cellfun(@(x)x.elapsed_s,fh)));
 fprintf('%d ',cellfun(@(x)x.end_visits,fh)); fprintf('\n');
end
% time-candidate enumeration health
tc=0; neg=0; pairs=0; refine=0; lcalls=0; prop=0; est=0; nat=0; fails=0;
for k=1:numel(S.construction_history)
 for j=1:numel(S.construction_history{k}.expansion_history)
  for t=1:numel(S.construction_history{k}.expansion_history{j}.time_search)
   r=S.construction_history{k}.expansion_history{j}.time_search{t}.screen;
   if ~isempty(r)&&isfield(r,'pheromone_queries')
    tc=tc+r.pheromone_queries; neg=neg+r.negative_guidance_hits; pairs=pairs+r.pairs_evaluated;
    refine=refine+r.time_refinements; lcalls=lcalls+r.lambert_calls; prop=prop+r.propagations;
    if isfield(r,'estimated_cost_rejections'), est=est+r.estimated_cost_rejections; end
    if isfield(r,'natural_windows'), nat=nat+r.natural_windows; end
    if isfield(r,'failures'), fails=fails+numel(r.failures); end
   end
  end
 end
end
fprintf('time candidates: queries %d | pairs %d | lambert %d | propagations %d | refinements %d\n', ...
 tc,pairs,lcalls,prop,refine);
fprintf('  estimated-cost rejections %d | natural windows %d | enumeration failures %d | negative hits %d\n', ...
 est,nat,fails,neg);
% joint calls
if ~isempty(S.history)
 rec=S.history;
 kinds=arrayfun(@(r)char(r.operator),cellfun(@(r)r,rec),'UniformOutput',false);
 [uk,~,ik]=unique(kinds); ck=accumarray(ik(:),1); [~,ok]=sort(ck,'descend');
 fprintf('joint/verify records %d\n',numel(rec));
 for k=reshape(ok,1,[]), fprintf('  %-40s %d\n',uk{k},ck(k)); end
 statuses=arrayfun(@(r)char(r.diagnostic.status),cellfun(@(r)r,rec),'UniformOutput',false);
 [us,~,is]=unique(statuses); cs=accumarray(is(:),1); [~,os]=sort(cs,'descend');
 for k=reshape(os,1,[]), fprintf('  status %-30s %d\n',us{k},cs(k)); end
end
% replans
if ~isempty(S.replan_history)
 rh=S.replan_history;
 fprintf('replans %d | total %.3f s | cost-driven %d\n',numel(rh),sum(cellfun(@(x)x.elapsed_s,rh)), ...
  sum(cellfun(@(x)isfield(x,'cost_driven')&&x.cost_driven,rh)));
 for k=1:numel(rh)
  r=rh{k}; np=0;
  if isfield(r,'attempts'), for a=1:numel(r.attempts), np=np+r.attempts{a}.proposals; end, end
  fprintf('  replan %d: status=%s cuts=[%s] proposals=%d elapsed=%.2f\n',k,r.status, ...
   num2str(r.cut_times_s(:).'),np,r.elapsed_s);
 end
end
% queues and rejections
fprintf('fresh_queue %d | recovery_pool %d | replan_queue %d | work_pool %d | refined_keys %d\n', ...
 numel(S.fresh_queue),numel(S.recovery_pool),numel(S.replan_queue),numel(S.work_pool),numel(S.refined_keys));
fprintf('policy_rejections %d | feedback events %d | prefix feedback %d\n', ...
 S.policy_rejection_count,numel(S.feedback_history),numel(S.prefix_feedback_history));
if ~isempty(S.best_partial)
 fprintf('best_prefix %d/35 J=%.6f t=%.1f\n',sum(S.best_partial.visited),S.best_partial.J,S.best_partial.t);
end
if isempty(S.best_rejected_complete), fprintf('no rejected complete candidate\n');
else, fprintf('best rejected complete %.9f km/s\n',S.best_rejected_complete.verification.total_dv_km_s); end
if isempty(S.elite), fprintf('no elite\n'); else, fprintf('ELITE %.9f km/s\n',S.elite.verification.total_dv_km_s); end
end
