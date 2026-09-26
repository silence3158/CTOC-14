# Natural-coast multi-flyby paper review

Date: 2026-09-26. Requested scope: assess the newly supplied paper's usefulness.
No algorithm edits, optimization batch, trajectory propagation experiment or
ATK execution. A small MATLAB calculation checked original target geometry.

## Source

An-yi Huang, Hong-xin Shen, Zhao Li, Cong Sun, Chao Sheng, Zheng-zhong Kuai,
*Global Optimization of Multi-Flyby Trajectories for Multi-Orbital-Plane
Constellations Inspection*, arXiv:2507.02943v1.
https://arxiv.org/abs/2507.02943

Local file: `literature/无机动自然滑行.pdf`.
SHA-256: B0FB930EBD395C3D512C26B9FC7D545D9B05E01367048316D772340EF13DC9AF.
Read Sections 2-5, algorithms, reported experiments and conclusion. Extracted
text: tmp/pdfs/natural_coast_20260926.txt. PDF page 10 rendered and inspected
to distinguish actual equation issues from extraction artifacts. No source
repository was identified in the supplied paper; no public implementation
was reviewed or imported. This paper is distinct from Zhang's epoch-pair DP.

## What it solves

The author assumes constellation planes populated by uniformly phased,
near-circular satellites. Instead of designing a separate transfer to every
satellite, select an inspection orbit with a small period offset. At repeated
near-perigee passages, different satellites arrive at the inspection location.
An entire plane can be inspected without maneuvers during the inspection arc.
Maneuvers remain necessary to enter the inspection orbit and transfer between
inspection orbits. The modeled mission starts on the first inspection orbit;
it does not include CTOC14's prescribed low initial-orbit departure cost.

For N uniformly spaced satellites with period T0, the positive-offset branch
uses inspection period T=(N+1)/N*T0. After each inspector revolution the
satellite constellation advances by an additional 2*pi/N modulo 2*pi; the
next appropriately phased satellite reaches the inspection location. Under
Kepler dynamics, the corresponding exact period-based semimajor axis is

    a = a0*(1+1/N)^(2/3).

The paper's Eq. (8), da approximately 2*a0/(3*N), is its first-order version.
This exact Kepler formula is explanatory, not a J2 encounter certificate.
Eccentricity places perigee near the target orbit; relative inclination/RAAN
and phase offsets control encounter geometry and J2 drift over the arc.
The interval from first to last of N visits is (N-1)*T, paper Eq. (7).

The three optimization layers are: (1) GA selects/orders orbital planes using
cheap transfer estimates, (2) DE reorders the selected planes and adjusts
arrival times and inspection-orbit offsets, (3) multi-impulse optimization
connects the fixed inspection orbits at the selected epochs. Algorithms are
not the main novelty; the physically constructed multi-target inspection
orbit is the useful reduction in search complexity.

## Actual assumptions and evidence boundary

| Item | Paper example | CTOC14 |
| --- | --- | --- |
| Targets | Structured LEO constellations; uniformly phased per plane | 35 heterogeneous targets |
| Dynamics | J2 averaged mean elements, Eq. (1) | Cartesian central gravity + J2 |
| Flyby distance | Less than 50 km | At most 1 km |
| Relative flyby speed | Less than 150 m/s | No speed-matching constraint |
| Duration | 90 days | At most 10 days |
| Objective | Maximize coverage under a 3 km/s budget | All 35 visits, minimize raw delta-V |
| First state | Starts on a selected inspection orbit | a0 in RE+[590,610] km, e0<0.001 |

Section 5.1 example: 22 satellites, approximately 210 km semimajor-axis offset,
1.457 days for the passive inspection arc, about 105 m/s relative speed. The
radial offset is deliberately 5 km, already outside CTOC14's flyby radius.
The validation uses Eq. (1), not our independent Cartesian J2 replay.

Section 5.2 reports 963 targets/32 planes, an estimated 3.673 km/s before
re-optimization and a final transfer-computed 2.98 km/s. The paper reports
CTOC13 first place with six spacecraft covering 5,516 targets/143 planes.
These are author-reported results under their conditions, not independently
reproduced CTOC14 outcomes. Their 25% exploratory budget relaxation is not
authorization to relax this project's 6.1 km/s search-pool threshold.

## Original target geometry check

MATLAB R2024a read data/ctoc14b_targets.csv. Plane normals were calculated from
the ORIGINAL Cartesian r/v, h=cross(r,v)/norm(cross(r,v)); pair angles use
acosd(clamp(dot(h1,h2))). This avoids treating near-equatorial RAAN differences
as physical plane angles. It does not propagate or prove flyby feasibility.

- Semimajor axes: 14441.238 to 42164.936 km.
- Eccentricities: 0.00031 to 0.69374.
- Inclinations: 0.12 to 63.48 degrees.
- Only 8 of 595 pairs have initial plane-angle difference below 1 degree.
- Only 3 pairs also have semimajor-axis difference below 100 km: (1,3),
  (1,2), (2,3), with plane angles 0.29687, 0.40049, 0.52911 degrees and
  semimajor-axis differences 0.458, 0.453, 0.911 km, respectively.

The thresholds 1 degree/100 km are descriptive diagnostics, not filters or
physical task restrictions. Target1/2/3 merit geometric examination as a
local candidate group, but are not an evenly phased coplanar constellation.
No group has been shown to admit a <=1 km natural multi-flyby arc.

Using the supplied a/e/i in the first-order secular estimate

    Omega_dot = -1.5*J2*(RE/(a*(1-e^2)))^2*sqrt(mu/a^3)*cos(i)

gives target rates roughly -0.570474 to -0.006966 degrees/day (mu=398600.4415,
RE=6378.1363, J2=1.08262668e-3). These are scale estimates using supplied
elements, not mean-element conversion or predictions replacing J2 integration.
For the supplied targets the estimated 10-day RAAN changes span about
0.07-5.70 degrees. Consequently 90-day LEO plane-alignment opportunities
cannot simply be assumed to recur within CTOC14's shorter mission. Inspector
drift can differ; this estimate is not an impossibility bound.

A first PowerShell scratch calculation of normals failed due to array
expression handling; all its pair-angle output was discarded. The numbers
above are from the corrected MATLAB calculation on Cartesian states.

## Formula caution

The supplied PDF itself has apparent notation/dimensional issues. Eq. (11)
prints sqrt(delta_v_flyby^2 - delta_v_t), with the tangential speed not
squared. A Euclidean two-component speed constraint requires its square.
Eq. (10)'s printed tangential expression adds terms of incompatible units.
Thus even where the geometric idea applies, formulas need to be rederived
and checked against relative-motion references and numerical dynamics; do
not blindly transcribe Algorithm 1. This observation does not invalidate
the author's reported numerical results.

## Transfer judgment

Useful: plan a natural arc covering a GROUP before optimizing transfers
between groups; choose period/phase to create timed encounters; keep small
plane offsets instead of unnecessarily matching each target's full orbit;
consider J2 drift when proposing windows. This gives a concrete literature
example for the user's 'maintain a basic state while accomplishing visits'
intuition, but only within appropriately structured inspection arcs.

Not established: that the full CTOC14 solution is stationkeeping, that its
35 targets can be partitioned into such analytical inspection planes, or
that this paper implies a sub-8/sub-6.1 km/s full trajectory.

In V3, v3Expand already includes zero-impulse coasts and random/target-guided
actions; v3GuidedChildren passes individual target-guided arcs to v3ApplyAction.
This supports discovering incidental extra visits but is not, by itself,
implementation of the paper's purpose-built resonant inspection orbit.
No exhaustive claim is made here about every V3 construction/refinement path.

A suitable later hypothesis test would construct one true natural J2 arc
with multiple encounter-time variables and no internal impulses, demanding
distance<=1 km for every proposed member and continuous altitude>=200 km.
Entering that arc and connecting its REAL exit velocity to the next segment
must be counted in full-task cost. A single successful short arc would prove
only a construction capability, not a complete cold-start mission.

The preceding Zhang delayed experiment retained one pulse per planned target.
Changing that wait time is different from deliberately removing internal
impulses to create a multi-target arc. A successful natural-arc constructor
could be represented as a multi-visit transition for outer search/DP, retaining
entry/exit time and state and the visited set; the exact two-epoch DP reduction
would generally no longer apply. This is an adaptation proposal, not an
implemented or validated mechanism.

Status: research archived; search model and code unchanged; no new run result.
