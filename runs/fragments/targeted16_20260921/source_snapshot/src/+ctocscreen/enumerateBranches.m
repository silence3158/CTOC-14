function [branches, info] = enumerateBranches(r1, r2, dtSeconds, mu, policy)
%ENUMERATEBRANCHES Enumerate two-body Lambert transfer branches.
%   [branches, info] = ctocscreen.enumerateBranches(r1, r2, dt, mu, policy)
%
% This offline solver brackets the universal-variable equation by pole intervals.
% It uses fzero and fminbnd for short-way and long-way
% geometries.  Positive-z roots are labelled by a revolution estimate; the
% estimate is a search classification, not a formal branch identifier.
%
% policy fields:
%   max_revolutions   non-negative integer, default 2
%   z_grid_points     legacy field, ignored by the bracketed solver
%   endpoint_tol_km   endpoint reconstruction tolerance, default 1e-3
%
% The solver returns an explicit no_branch status when nothing is found.
% It never silently drops failed roots.

if nargin < 5 || isempty(policy)
    policy = struct();
end
policy = fillDefault(policy, 'max_revolutions', 2);
policy = fillDefault(policy, 'z_grid_points', 800);
policy = fillDefault(policy, 'endpoint_tol_km', 1e-3);
validateattributes(policy.max_revolutions,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(policy.endpoint_tol_km,{'numeric'},{'scalar','positive','finite'});
validateattributes(mu,{'numeric'},{'scalar','positive','finite'});
validateattributes(r1,{'numeric'},{'real','finite','numel',3});
validateattributes(r2,{'numeric'},{'real','finite','numel',3});

r1 = r1(:).';
r2 = r2(:).';
if numel(r1) ~= 3 || numel(r2) ~= 3
    error('ctocscreen:enumerateBranches:badVectors', ...
        'r1 and r2 must be 3-vectors.');
end
if ~(isscalar(dtSeconds) && isfinite(dtSeconds) && dtSeconds > 0)
    error('ctocscreen:enumerateBranches:badDt', ...
        'dtSeconds must be a positive finite scalar.');
end

r1n = norm(r1);
r2n = norm(r2);
cosd = dot(r1, r2) / (r1n*r2n);
cosd = max(-1, min(1, cosd));
base = acos(cosd);

branches = emptyBranchStruct();
info = struct();
info.status = 'ok';
info.message = '';
info.r1n_km = r1n;
info.r2n_km = r2n;
info.transfer_angle_short_rad = base;
info.policy = policy;

if abs(sin(base)) < 1e-10
    info.status = 'collinear_degenerate';
    info.message = 'Transfer geometry is numerically collinear.';
    return;
end


geometry = struct( ...
    'name', {'short', 'long'}, ...
    'dnu', {base, 2*pi - base});

allBranches = emptyBranchStruct();
nRoots = 0;
for g = 1:numel(geometry)
    A = sin(geometry(g).dnu) * sqrt(r1n*r2n/(1 - cosd));
    zroots = bracketRoots(r1n,r2n,A,mu,dtSeconds,policy.max_revolutions);
    nRoots = nRoots + numel(zroots);
    for j = 1:numel(zroots)
        [branch, ok] = makeBranch(zroots(j), geometry(g).name, ...
            geometry(g).dnu, A, r1, r2, r1n, r2n, dtSeconds, mu, policy);
        if ok
            rev=branch.rev_estimate;
            earlier=zroots(1:j-1);
            side=sum(floor(sqrt(max(0,earlier))/(2*pi))==rev);
            branch.branch_id=1+4*rev+2*(g-1)+side;
            allBranches(end+1) = branch; %#ok<AGROW>
        end
    end
end

info.n_roots = nRoots;
info.n_branches = numel(allBranches);
info.n_rejected_roots = nRoots-numel(allBranches);
if isempty(allBranches)
    info.status = 'no_branch';
    info.message = 'No valid Lambert branch found on the scanned z intervals.';
    branches = allBranches;
    return;
end

% Stable ordering: revolution estimate, |z|, geometry name.
key = [ [allBranches.rev_estimate].', abs([allBranches.z]).', ...
    double([allBranches.geometry].' == "short") ];
[~, order] = sortrows(key, [1 2 3]);
allBranches = allBranches(order);

branches = allBranches;
info.status = 'ok';
end

function roots = bracketRoots(r1,r2,A,mu,dt,maxrev)
% Each positive pole interval has one minimum and at most two roots.
fun=@(z) timeValue(z,r1,r2,A,mu)-dt;
opt=optimset('Display','off','TolX',2e-12);
roots=[];
hi=(2*pi)^2-1e-5; lo=0;
while fun(lo)>0 && lo>-1e5, lo=2*lo-1; end
if fun(lo)<=0 && fun(hi)>0
 roots(end+1)=fzero(fun,[lo hi],opt);
end
for m=1:maxrev
 lo=(2*pi*m)^2+1e-5; hi=(2*pi*(m+1))^2-1e-5;
 [zm,fm]=fminbnd(fun,lo,hi,opt);
 if fm < 0
  roots(end+1)=fzero(fun,[lo zm],opt); %#ok<AGROW>
  roots(end+1)=fzero(fun,[zm hi],opt); %#ok<AGROW>
 elseif abs(fm)<1e-8
  roots(end+1)=zm; %#ok<AGROW>
 end
end
end

function t=timeValue(z,r1,r2,A,mu)
[C,S]=ctocscreen.stumpff(z);
y=r1+r2+A*(z*S-1)/sqrt(C);
if y<=0, t=0; return; end
t=((y/C)^1.5*S+A*sqrt(y))/sqrt(mu);
end

% -------------------------------------------------------------------------
function b = branchPrototype()
b = struct( ...
    'branch_id', 0, ...
    'geometry', "", ...
    'transfer_angle_rad', NaN, ...
    'rev_estimate', 0, ...
    'z', NaN, ...
    'A', NaN, ...
    'y', NaN, ...
    'v_depart', nan(1,3), ...
    'v_arrive', nan(1,3), ...
    'end_residual_km', NaN, ...
    'status', "");
end

function b = emptyBranchStruct()
b = repmat(branchPrototype(), 0, 1);
end

function policy = fillDefault(policy, field, value)
if ~isfield(policy, field) || isempty(policy.(field))
    policy.(field) = value;
end
end

function [branch, ok] = makeBranch(z, geomName, dnu, A, r1, r2, r1n, r2n, dt, mu, policy)
ok = false;
branch = branchPrototype();
[C, S] = ctocscreen.stumpff(z);
if ~isfinite(C) || C <= 0
    return;
end
y = r1n + r2n + A*(z*S - 1)/sqrt(C);
if ~isfinite(y) || y <= 0
    return;
end
f = 1 - y/r1n;
g = A*sqrt(y/mu);
gdot = 1 - y/r2n;
if ~isfinite(g) || abs(g) < 1e-14
    return;
end
vdep = (r2 - f*r1) / g;
varr = (gdot*r2 - r1) / g;
if any(~isfinite(vdep)) || any(~isfinite(varr))
    return;
end
[stateChk, infChk] = ctocscreen.propagateTwoBody([r1, vdep], dt, mu);
if ~strcmp(infChk.status, 'ok')
    return;
end
res = norm(stateChk(1:3) - r2);
if ~isfinite(res) || res > policy.endpoint_tol_km
    return;
end
rev = 0;
if z > 0
    rev = floor(sqrt(z)/(2*pi));
end
branch.branch_id = 0;
branch.geometry = string(geomName);
branch.transfer_angle_rad = dnu;
branch.rev_estimate = rev;
branch.z = z;
branch.A = A;
branch.y = y;
branch.v_depart = vdep;
branch.v_arrive = varr;
branch.end_residual_km = res;
branch.status = "ok";
ok = true;
end
