function batch=run_batch_screening(mode,workers,runId,resume)
%RUN_BATCH_SCREENING Example: run_batch_screening('batch',4,'night01',false)
if nargin<1, mode='smoke'; end
if nargin<2, workers=0; end
if nargin<4, resume=false; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
config=ctocscreen.defaultConfig(mode); config.workers=workers;
if nargin>=3 && strlength(string(runId))>0, config.run_id=string(runId); end
problem=ctocscreen.loadProblem();
batch=ctocscreen.runBatch(problem,config,resume);
end
