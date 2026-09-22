function problem = loadProblem(csvPath, sourceMetadata)
%LOADPROBLEM Load and validate CTOC14B target data for screening.
%   problem = ctocscreen.loadProblem(csvPath)
%   problem = ctocscreen.loadProblem(csvPath, sourceMetadata)
%
% The CSV is data/ctoc14b_targets.csv.  Positions are km, velocity
% components are km/s, angles are degrees.  Internal units are km, s, rad,
% km/s.  This function performs format-level validation only; coordinate
% frame consistency with the original ATK scenario remains a separate task.

isEmptyPath = (nargin < 1) || isempty(csvPath);
if ~isEmptyPath && (ischar(csvPath) || isstring(csvPath))
    isEmptyPath = strlength(string(csvPath)) == 0;
end
if isEmptyPath
    thisFile = mfilename('fullpath');
    pkgDir = fileparts(thisFile);
    srcDir = fileparts(pkgDir);
    simDir = fileparts(srcDir);
    csvPath = fullfile(simDir, 'data', 'ctoc14b_targets.csv');
end
if nargin < 2 || isempty(sourceMetadata)
    sourceMetadata = struct();
end
if ~isfile(csvPath)
    error('ctocscreen:loadProblem:missingFile', ...
        'Target CSV not found: %s', csvPath);
end

T = readtable(csvPath, 'VariableNamingRule', 'preserve', 'TextType', 'string');
T=sortrows(T,1);
if width(T) < 25
    error('ctocscreen:loadProblem:badColumns', ...
        'Target CSV must have at least 25 columns, found %d.', width(T));
end

ids = double(T{:,1});
names = string(T{:,2});
epoch = string(T{:,3});
states = [ ...
    double(T{:,18}), double(T{:,19}), double(T{:,20}), ...
    double(T{:,21}), double(T{:,22}), double(T{:,23})];
rnorm_csv = double(T{:,24});
vnorm_csv = double(T{:,25});
n = height(T);

if n ~= 35
    error('ctocscreen:loadProblem:targetCount', ...
        'Expected 35 targets, found %d.', n);
end
if numel(unique(ids)) ~= n || any(sort(ids) ~= (1:n).')
    error('ctocscreen:loadProblem:targetIds', ...
        'Target IDs must be unique and equal to 1..35.');
end
expectedNames = "Target" + string((1:n).');
if any(names ~= expectedNames)
    error('ctocscreen:loadProblem:targetNames', ...
        'Target names must be Target1..Target35.');
end
if numel(unique(epoch)) ~= 1 || epoch(1)~="2035-01-01 12:00:00"
    error('ctocscreen:loadProblem:epoch', ...
        'All targets must share one epoch.');
end
if any(~isfinite(states), 'all')
    error('ctocscreen:loadProblem:finiteState', ...
        'Target Cartesian states must be finite.');
end
if any(vecnorm(states(:,1:3), 2, 2) <= 0)
    error('ctocscreen:loadProblem:zeroPosition', ...
        'Target position vectors must be non-zero.');
end
if any(vecnorm(states(:,4:6), 2, 2) <= 0)
    error('ctocscreen:loadProblem:zeroVelocity', ...
        'Target velocity vectors must be non-zero.');
end

rnorm = vecnorm(states(:,1:3), 2, 2);
vnorm = vecnorm(states(:,4:6), 2, 2);
if max(abs(rnorm - rnorm_csv)) > 1e-6
    warning('ctocscreen:loadProblem:rnormMismatch', ...
        '|r| column differs from Cartesian position norm.');
end
if max(abs(vnorm*1000 - vnorm_csv)) > 1e-3
    warning('ctocscreen:loadProblem:vnormMismatch', ...
        '|v| column differs from Cartesian velocity norm times 1000.');
end

problem = struct();
problem.target_ids = ids(:).';
problem.target_names = names(:).';
problem.epoch_utc = epoch(1);
problem.frame = string(getfield_default(sourceMetadata, 'frame', "unknown"));
problem.states0 = states;
problem.mu_km3_s2 = 398600.4415;
problem.re_km = 6378.1363;
problem.horizon_s = 864000;
problem.source_path = string(csvPath);
problem.source_metadata = sourceMetadata;
problem.input_hash = hashFile(csvPath);
problem.model_id = "two_body_screen_v2";
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
scenario=fullfile(root,'problem','ctoc14乙题','CTOC14B.atk');
gravity=fullfile(root,'problem','ctoc14乙题','ATK-CTOC14','AstroData','Earth','EGM96.grv');
problem.source_audit=ctocscreen.auditSource(problem,scenario,gravity);
problem.frame="ATK_native_inertial";
end

function v = getfield_default(s, f, d)
if isfield(s, f)
    v = s.(f);
else
    v = d;
end
end

function h = hashFile(fname)
fid = fopen(fname, 'r');
if fid < 0
    error('ctocscreen:loadProblem:hashOpen', ...
        'Cannot open file for hashing: %s', fname);
end
bytes = fread(fid, Inf, '*uint8');
fclose(fid);
md = java.security.MessageDigest.getInstance('SHA-256');
digest = typecast(md.digest(typecast(bytes,'int8')), 'uint8');
h = lower(reshape(dec2hex(digest, 2).', 1, []));
end

