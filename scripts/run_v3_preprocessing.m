function folder = run_v3_preprocessing(label,options)
%RUN_V3_PREPROCESSING Build nominal J2 cache; does not invoke ATK.
if nargin<1, label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
if nargin<2, options=struct(); end
assert(~isempty(regexp(label,'^[A-Za-z0-9_-]+$','once')), ...
 'ctocscreen:v3:label','Label must be a simple unique directory name.');
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
folder=fullfile(sim,'runs','v3','preprocessing',label);
ctocscreen.v3BuildTargetEphemeris(folder,options);
fprintf('Cache and report: %s\n',folder);
end
