function state = loadCheckpoint(filePath)
%LOADCHECKPOINT Load a MAT checkpoint if it exists.
%   state = ctocscreen.loadCheckpoint(filePath)

if ~isfile(filePath)
    error('ctocscreen:checkpoint:missing','Checkpoint file does not exist.');
end
S = load(filePath);
required={'signature','iter','archive','history','streamState','config','elapsed_s','refinements'};
assert(all(isfield(S,required)),'ctocscreen:checkpoint:format', ...
 'Checkpoint is incomplete or from an incompatible pre-v2 run.');
state = S;
end
