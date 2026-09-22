function saveCheckpoint(filePath, state)
%SAVECHECKPOINT Atomically write a MAT checkpoint.
%   ctocscreen.saveCheckpoint(filePath, state)

if nargin < 2
    error('ctocscreen:saveCheckpoint:missingState', 'state is required.');
end
parentDir = fileparts(filePath);
if ~isempty(parentDir) && ~isfolder(parentDir)
    mkdir(parentDir);
end
tmpPath = [filePath, '.tmp'];
save(tmpPath, '-struct', 'state', '-v7');
% Preserve the previous valid generation even if replacement is interrupted.
if isfile(filePath), copyfile(filePath,[filePath '.bak'],'f'); end
[ok,msg]=movefile(tmpPath, filePath, 'f');
if ~ok, error('ctocscreen:checkpoint:replace','%s',msg); end
end
