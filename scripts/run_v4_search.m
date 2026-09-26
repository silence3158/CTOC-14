function folder=run_v4_search(label,budget_s,seed,options,ephemerisFile)
%RUN_V4_SEARCH Empty-history A+B cold search with independent final verification.
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
if nargin<2, budget_s=300; end
if nargin<3, seed=888; end
if nargin<4, options=struct(); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
if nargin<5, ephemerisFile=fullfile(sim,'runs','v3','preprocessing','targets_20260923_release','target_ephemeris.mat'); end
assert(~isempty(regexp(label,'^[A-Za-z0-9_-]+$','once')),'ctocscreen:v4:label','Invalid run label.');
assert(license('test','Optimization_Toolbox'),'ctocscreen:v4:license','Optimization Toolbox is required.');
options.budget_s=budget_s; options.seed=seed; c=ctocscreen.v4.defaults(options);
eph=ctocscreen.v3LoadTargetEphemeris(ephemerisFile);
assert(~isfield(eph,'synthetic')||~eph.synthetic,'ctocscreen:v4:dataset','Competition run cannot use synthetic data.');
folder=fullfile(sim,'runs','v4','search',label); ctocscreen.v4.search(folder,eph,c);
fprintf('V4 report: %s\n',fullfile(folder,'report.md'));
end
