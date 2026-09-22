% Compatibility entry point; small serial smoke run.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
problem=ctocscreen.loadProblem();
config=ctocscreen.defaultConfig('smoke');
result=ctocscreen.runSearch(problem,config);
disp(result.summary);