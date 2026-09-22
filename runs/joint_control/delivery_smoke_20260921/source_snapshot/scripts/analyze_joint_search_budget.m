function report = analyze_joint_search_budget()
%ANALYZE_JOINT_SEARCH_BUDGET Read-only trajectory diagnostics, no optimization.
root = fileparts(fileparts(mfilename('fullpath')));
source = fullfile(root,'runs','fragments','diverse16b_20260921','elite.mat');
d = load(source,'elite','p');
s = d.elite.schedule;
dv = vecnorm(s.delta_v_km_s,2,2);
[sorted, order] = sort(dv,'descend');
burns = table((1:numel(dv))',s.maneuver_times_s(:)/86400,dv, ...
    'VariableNames',{'burn','time_day','dv_km_s'});
x = d.p.states0;
mu = d.p.mu_km3_s2;
a = 1./(2./vecnorm(x(:,1:3),2,2)-sum(x(:,4:6).^2,2)/mu);
period = 2*pi*sqrt(a.^3/mu)/3600;
targets = table((1:size(x,1))',a,period, ...
    'VariableNames',{'target','semimajor_axis_km','period_h'});
report = struct('source',source,'total_dv_km_s',sum(dv), ...
    'first_dv_km_s',dv(1),'top5_dv_km_s',sum(sorted(1:5)), ...
    'top10_dv_km_s',sum(sorted(1:10)), ...
    'below_50m_s_count',sum(dv<.05),'below_50m_s_sum_km_s',sum(dv(dv<.05)), ...
    'burn_rank',order,'burns',burns,'targets',targets, ...
    'scope','Diagnostics of archived two-body elite; no new feasibility verification');
out = fullfile(root,'runs','analysis','joint_search_20260921');
if ~isfolder(out), mkdir(out); end
save(fullfile(out,'budget.mat'),'report');
writetable(burns,fullfile(out,'burn_budget.csv'));
writetable(targets,fullfile(out,'target_periods.csv'));
disp(report);
end
