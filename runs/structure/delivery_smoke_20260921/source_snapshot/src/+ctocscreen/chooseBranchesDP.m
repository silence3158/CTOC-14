function path = chooseBranchesDP(branchSets, vStart)
%CHOOSEBRANCHESDP Dynamic program over one branch per leg.
%   path = ctocscreen.chooseBranchesDP(branchSets, vStart)
%
% branchSets{k} is a struct array returned by enumerateBranches for leg k.
% The DP minimizes sum(norm(v_depart[k] - v_arrive[k-1])) with
% v_arrive[0] = vStart.  It is exact for the nominal fixed-boundary branch
% sets, not for the full continuous problem.

N = numel(branchSets);
if N < 1
    error('ctocscreen:chooseBranchesDP:empty', 'branchSets is empty.');
end
maxB = max(1,max(cellfun(@numel, branchSets)));
D = inf(N, maxB);
prev = zeros(N, maxB);

for b = 1:numel(branchSets{1})
    D(1,b) = norm(branchSets{1}(b).v_depart - vStart);
end

for k = 2:N
    for b = 1:numel(branchSets{k})
        vb = branchSets{k}(b).v_depart;
        bestCost = inf;
        bestPrev = 0;
        for a = 1:numel(branchSets{k-1})
            if ~isfinite(D(k-1,a))
                continue;
            end
            c = D(k-1,a) + norm(vb - branchSets{k-1}(a).v_arrive);
            if c < bestCost
                bestCost = c;
                bestPrev = a;
            end
        end
        D(k,b) = bestCost;
        prev(k,b) = bestPrev;
    end
end

[bestCost, bestB] = min(D(N,:));
path = struct();
path.status = 'ok';
path.nominal_dv_km_s = bestCost;
path.branch_indices = zeros(1,N);
path.branch_ids = zeros(1,N);
path.selected = cell(1,N);
path.message = "";
path.failure_leg = 0;

if ~isfinite(bestCost)
    path.status = 'no_branch';
    path.message = 'No complete branch chain exists.';
    path.failure_leg=find(cellfun(@isempty,branchSets),1);
    return;
end

idx = zeros(1,N);
idx(N) = bestB;
for k = N:-1:2
    idx(k-1) = prev(k, idx(k));
    if idx(k-1) == 0
        path.status = 'no_branch';
        path.message = 'Backtracking failed.';
        return;
    end
end
path.branch_indices = idx;
for k = 1:N
    path.branch_ids(k) = branchSets{k}(idx(k)).branch_id;
    path.selected{k} = branchSets{k}(idx(k));
end
end

