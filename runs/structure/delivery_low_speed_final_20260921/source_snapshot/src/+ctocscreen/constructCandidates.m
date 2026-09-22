function candidates = constructCandidates(problem, config, stream, count)
%CONSTRUCTCANDIDATES Generate randomized initial screening candidates.
%   candidates = ctocscreen.constructCandidates(problem, config, stream, count)

if nargin < 3 || isempty(stream)
    stream = RandStream('mt19937ar', 'Seed', 1);
end
if nargin < 4 || isempty(count)
    count = 1;
end
if ~isfield(config, 'construct') || isempty(config.construct)
    config.construct = struct();
end
c = config.construct;
if ~isfield(c, 'wait_max_s'), c.wait_max_s = 3600; end
if ~isfield(c, 'tof_total_fraction'), c.tof_total_fraction = [0.2, 0.8]; end
if ~isfield(c, 'min_tof_s'), c.min_tof_s = 1.0; end

N = 35;
re = problem.re_km;
aMin = re + 590;
aMax = re + 610;
candidates = [];

for i = 1:count
    a0 = aMin + (aMax - aMin)*rand(stream);
    e0 = 0.0009 * rand(stream);
    omega0 = 2*pi*rand(stream);
    ec0 = e0 * cos(omega0);
    es0 = e0 * sin(omega0);
    cosi = 1 - 2*rand(stream);
    i0 = acos(max(-1, min(1, cosi)));
    Omega0 = 2*pi*rand(stream);
    u0 = 2*pi*rand(stream);
    q0 = [a0, ec0, es0, i0, Omega0, u0];

    order = randperm(stream, N);
    wait_s = c.wait_max_s * rand(stream);
    window = problem.horizon_s - wait_s;
    frac = c.tof_total_fraction(1) ...
        + (c.tof_total_fraction(2) - c.tof_total_fraction(1))*rand(stream);
    totalTof = max(N*c.min_tof_s, window*frac);
    if totalTof > window
        totalTof = window;
    end
    base = 0.2 + rand(stream, 1, N);
    base = base / sum(base);
    tof_s = c.min_tof_s + (totalTof - N*c.min_tof_s)*base;

    cand = ctocscreen.makeCandidate(q0, order, wait_s, tof_s, ...
        zeros(1,N), randi(stream, 2^31-1), "", "random_construct");
    if i == 1
        candidates = cand;
    else
        candidates(i) = cand; %#ok<AGROW>
    end
end
end

