function folder=run_v3_search(label,budget_s,seed,options,ephemerisFile)
%RUN_V3_SEARCH Complete V3 nominal-J2 search, independent of old archives.
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
if nargin<2, budget_s=60; end
if nargin<3, seed=1; end
if nargin<4, options=struct(); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
if nargin<5, ephemerisFile=fullfile(sim,'runs','v3','preprocessing','targets_20260923_release','target_ephemeris.mat'); end
assert(~isempty(regexp(label,'^[A-Za-z0-9_-]+$','once')),'ctocscreen:v3:label','Invalid run label.');
options.budget_s=budget_s; options.seed=seed; c=ctocscreen.v3Defaults(options);
eph=ctocscreen.v3LoadTargetEphemeris(ephemerisFile);
folder=fullfile(sim,'runs','v3','search',label);
ctocscreen.v3Search(folder,eph,c);
fprintf('V3 run: %s\n',folder);
end
