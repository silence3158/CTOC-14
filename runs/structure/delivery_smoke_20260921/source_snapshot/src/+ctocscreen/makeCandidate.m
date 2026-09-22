function candidate = makeCandidate(initial_q, order, wait_s, tof_s, branch_ids, seed, parent_id, generation_method)
%MAKECANDIDATE Construct a screening candidate structure.
%   candidate = ctocscreen.makeCandidate(initial_q, order, wait_s, tof_s)
%   candidate = ctocscreen.makeCandidate(..., branch_ids, seed, parent_id, method)

if nargin < 5 || isempty(branch_ids)
    branch_ids = zeros(1, numel(order));
end
if nargin < 6 || isempty(seed)
    seed = NaN;
end
if nargin < 7 || isempty(parent_id)
    parent_id = "";
end
if nargin < 8 || isempty(generation_method)
    generation_method = "manual";
end

initial_q = initial_q(:).';
order = order(:).';
wait_s = double(wait_s);
tof_s = tof_s(:).';
branch_ids = branch_ids(:).';

if numel(initial_q) ~= 6
    error('ctocscreen:makeCandidate:initialQ', ...
        'initial_q must have 6 elements.');
end
N = numel(order);
if N < 1 || numel(tof_s) ~= N || numel(branch_ids) ~= N
    error('ctocscreen:makeCandidate:sizeMismatch', ...
        'order, tof_s and branch_ids must have the same length.');
end
if any(sort(order) ~= (1:N))
    error('ctocscreen:makeCandidate:order', ...
        'order must be a permutation of 1..N.');
end
if ~(isscalar(wait_s) && isfinite(wait_s) && wait_s >= 0)
    error('ctocscreen:makeCandidate:wait', ...
        'wait_s must be a non-negative finite scalar.');
end
if any(~isfinite(tof_s)) || any(tof_s <= 0)
    error('ctocscreen:makeCandidate:tof', ...
        'tof_s must contain positive finite values.');
end

candidate = struct();
candidate.initial_q = initial_q;
candidate.order = order;
candidate.wait_s = wait_s;
candidate.tof_s = tof_s;
candidate.branch_ids = branch_ids;
candidate.seed = seed;
candidate.parent_id = string(parent_id);
candidate.generation_method = string(generation_method);
candidate.encoding_version = 1;
candidate.candidate_id = "";
end
