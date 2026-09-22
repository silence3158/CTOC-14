function archive = updateArchive(archive, candidate, evaluation, maxSize)
%UPDATEARCHIVE Keep the best complete two-body-screened candidates.
%   This first-stage archive is ordered by total_dv_km_s.  It is not a
%   diversity archive; that is a later search-stage feature.

if nargin < 4 || isempty(maxSize)
    maxSize = 16;
end
if isempty(archive)
    archive = repmat(struct('candidate', [], 'evaluation', [], ...
        'total_dv_km_s', inf), 0, 1);
end
if ~strcmp(evaluation.status, 'two_body_verified') || ...
        evaluation.reached_planned_count~=35 || evaluation.unique_visit_count~=35
    return;
end
if ~isfinite(evaluation.total_dv_km_s)
    return;
end

entry = struct();
entry.candidate = candidate;
entry.evaluation = evaluation;
entry.total_dv_km_s = evaluation.total_dv_km_s;
for k=1:numel(archive)
    old=archive(k).candidate;
    if isequal(old.initial_q,candidate.initial_q) && isequal(old.order,candidate.order) && ...
       isequal(old.tof_s,candidate.tof_s) && old.wait_s==candidate.wait_s && ...
       isequal(old.branch_ids,candidate.branch_ids)
        return;
    end
end
archive(end+1) = entry;
[~, order] = sort([archive.total_dv_km_s], 'ascend');
archive = archive(order);
if numel(archive) > maxSize
    archive = archive(1:maxSize);
end
end
