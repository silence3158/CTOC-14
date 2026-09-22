function [state, info] = initialState(q0, mu, re)
%INITIALSTATE Convert the six initial-orbit variables to a Cartesian state.
%   q0 = [a0, ec0, es0, i0, Omega0, u0]
%   a0 in km, ec0/es0 dimensionless, angles in rad.
%   u0 = omega0 + nu0 is the argument of latitude at epoch.
%
%   [state, info] = ctocscreen.initialState(q0, mu, re)
%   state = [x y z vx vy vz], km and km/s.

if nargin < 2 || isempty(mu)
    mu = 398600.4415;
end
if nargin < 3 || isempty(re)
    re = 6378.1363;
end
q0 = q0(:).';
if numel(q0) ~= 6 || any(~isfinite(q0))
    error('ctocscreen:initialState:badQ', ...
        'q0 must be a finite 6-vector.');
end

a = q0(1);
ec = q0(2);
es = q0(3);
inc = q0(4);
Omega = q0(5);
u = q0(6);

if ~(a > 0)
    error('ctocscreen:initialState:semiMajor', ...
        'Semi-major axis must be positive.');
end

e = hypot(ec, es);
if e < 1e-12
    omega = 0;
    nu = mod(u, 2*pi);
else
    omega = atan2(es, ec);
    nu = mod(u - omega, 2*pi);
end

if e >= 1
    error('ctocscreen:initialState:ellipticOnly', ...
        'Initial orbit must be elliptic (e < 1).');
end
p = a * (1 - e^2);
if ~(p > 0)
    error('ctocscreen:initialState:badP', ...
        'Semi-latus rectum must be positive.');
end

den = 1 + e*cos(nu);
if den <= 0
    error('ctocscreen:initialState:badRadius', ...
        'Invalid true anomaly for elliptic orbit.');
end
r = p / den;
r_pf = [r*cos(nu); r*sin(nu); 0];
v_pf = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];

cO = cos(Omega); sO = sin(Omega);
ci = cos(inc);   si = sin(inc);
cw = cos(omega); sw = sin(omega);
R = [ ...
    cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci,  sO*si; ...
    sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si; ...
    sw*si,             cw*si,             ci];

r_eci = R * r_pf;
v_eci = R * v_pf;
state = [r_eci.', v_eci.'];

info = struct();
info.e = e;
info.omega = omega;
info.nu = nu;
info.p = p;
info.re = re;
end
