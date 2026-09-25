function report = report_beam_v1(label, baselineDv, referenceDv)
%REPORT_BEAM_V1 V1-style report for a beam batch: only the figures/table that the
%V1 scripts produced (tmp/archive_expanded.m and tmp/archive_campaign16.m).
%
%   convergence.png   x = outer iteration, y = total delta-V, one line per task,
%                     dashed previous-baseline line, dotted user reference line
%                     (same style as archive_expanded.m)
%   wave_results.png  x = task index, y = best INDEPENDENTLY VERIFIED delta-V,
%                     with the user reference line (same style as
%                     archive_campaign16.m's wave plot; x is the task here
%                     because a beam batch is organised as tasks, not waves)
%   trajectory.csv    leg,target,departure_s,arrival_s,dvx,dvy,dvz,dv,independent_error
%                     for the best independently verified solution
%
%   report = report_beam_v1('beam_v1_a')
%   report = report_beam_v1('beam_v1_a', 16.287951653319, 7)
%
% Only reads runs/screening/<label>/batch.mat. Never touches src.

if nargin < 1 || isempty(label)
    error('report_beam_v1:label','Provide the run label, e.g. report_beam_v1(''beam_v1_a'').');
end
if nargin < 2 || isempty(baselineDv)
    baselineDv = 16.287951653319;   % V1 campaign16_20260921 verified baseline (km/s)
end
if nargin < 3 || isempty(referenceDv)
    referenceDv = 7;                % user-reported (not independently verified) reference
end
root = fileparts(fileparts(mfilename('fullpath')));
cd(root); addpath('src'); addpath('beamsearch');
folder = fullfile('runs','screening',label);
assert(isfolder(folder),'report_beam_v1:missing','No such run folder: %s',folder);
S = load(fullfile(folder,'batch.mat'));
batch = S.batch; p = batch.problem;
n = numel(batch.results);

% ---- convergence.png (style of tmp/archive_expanded.m) ---------------------
fig = figure('Visible','off','Position',[100 100 1000 550]); hold on;
for j = 1:n
    o = batch.results{j}; h = o.history;
    if isempty(h), continue; end
    plot([h.iter],[h.best_dv],'LineWidth',1.5, ...
        'DisplayName',sprintf('Task %d (seed %d)', j, o.config.master_seed));
end
xlabel('Outer iteration'); ylabel('Total delta-V (km/s)'); grid on;
yline(baselineDv,'--','Previous verified baseline','HandleVisibility','off');
yline(referenceDv,':','User-reported reference','HandleVisibility','off');
legend('Location','northeast');
title(sprintf('Beam V1 search (width %d, %d tasks): screened objective histories', ...
    batch.config.beam_width, n));
exportgraphics(fig,fullfile(folder,'convergence.png'),'Resolution',150); close(fig);

% ---- wave_results.png (style of tmp/archive_campaign16.m) ------------------
verified = nan(1,n); screened = nan(1,n); seeds = zeros(1,n);
for j = 1:n
    o = batch.results{j};
    seeds(j) = o.config.master_seed;
    screened(j) = o.summary.best_screened_dv_km_s;
    if isfinite(o.summary.best_verified_dv_km_s)
        verified(j) = o.summary.best_verified_dv_km_s;
    end
end
fig = figure('Visible','off','Position',[100 100 900 450]); hold on;
plot(1:n,verified,'-o','LineWidth',1.5,'DisplayName','independently verified');
plot(1:n,screened,'-s','LineWidth',1.0,'DisplayName','screened (not verified)');
yline(baselineDv,'--','V1 baseline','HandleVisibility','off');
yline(referenceDv,':','User-reported reference: 7','HandleVisibility','off');
grid on; xlabel('Task'); ylabel('Best total delta-V (km/s)');
title('Beam V1 batch: screened and verified results by task');
xticks(1:n); legend('Location','best');
exportgraphics(fig,fullfile(folder,'wave_results.png'),'Resolution',150); close(fig);

% ---- trajectory.csv of the best verified solution --------------------------
report = struct('folder',folder,'baseline_dv_km_s',baselineDv, ...
    'reference_dv_km_s',referenceDv,'best',[],'seeds',seeds);
if isempty(batch.export_archive)
    fprintf('No independently verified candidate in this batch; only convergence.png written.\n');
    return;
end
e = batch.export_archive(1); e = e(1);
report.best = e;
c = e.candidate; ev = e.evaluation; a = ev.independent_report;
leg = (1:numel(c.order))';
target = c.order';
departure_s = ev.depart_times_s';
arrival_s = ev.arrive_times_s';
dvx_km_s = ev.delta_v(:,1); dvy_km_s = ev.delta_v(:,2); dvz_km_s = ev.delta_v(:,3);
dv_km_s = vecnorm(ev.delta_v,2,2);
independent_error_km = a.endpoint_errors_km';
writetable(table(leg,target,departure_s,arrival_s,dvx_km_s,dvy_km_s,dvz_km_s, ...
    dv_km_s,independent_error_km), fullfile(folder,'trajectory.csv'));

fprintf(['FINAL_BEAM best %.12f gain_vs_V1_baseline %.6f%% visits %d ' ...
    'max_error_m %.9f min_alt %.9f days %.9f burns %d\n'], ...
    a.total_dv_km_s, 100*(1-a.total_dv_km_s/baselineDv), a.unique_visit_count, ...
    1000*a.max_endpoint_error_km, a.min_altitude_km, ev.duration_s/86400, numel(c.order));
end
