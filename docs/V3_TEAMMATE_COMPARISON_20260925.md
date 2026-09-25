# Teammate strategy comparison, 2026-09-25

## Scope and evidence boundary

The user supplied `tail27_polished.atk` for diagnostic comparison and confirmed
that the reported result is raw delta-V slightly below 10 km/s, with F2 about
9.72. That reported result has NOT been reproduced here. No ATK was launched,
no optimization search was run, and no teammate or historical schedule was
provided to the cold-start search.

Input SHA-256:
`2D4506B263E710291F814A9DC5193A47482F5A9C89D1BBB2E20DE9438B5A8532`.

The XML is a strategy specification with inconsistent cached states, not a
usable fixed-impulse record. A local reconstruction is a DIFFERENT result until
ATK branch semantics and actual executed impulses have been reconciled.

## XML findings

- One initial-state segment, 30 coast segments, 35 single-departure-impulse
  Lambert target segments, all active; all 35 target IDs occur once.
- Lambert segments specify perturbations enabled, zero revolutions,
  `MinorArc=-1`, `ShortWay=0`. Exact correspondence to our branch selection has
  not been established.
- Configured initial orbit: a=6987.1363 km (609 km above reference radius),
  essentially circular, inclination 45 degrees, RAAN 30 degrees.
- Summed duration: 828759.70572 s, or 9.592126 days. The sum agrees with configured
  Lambert arrival UTCs to rounding precision.
- Total explicit coasting: 1.695481 days. Twenty waits exceed 60 s, including the
  initial wait; nineteen such waits occur after a visit. Longest wait: 8.690620 h.
- Median Lambert flight duration: 5.724533 h.
- Raw target initial-state components match the local input to at most
  4.95e-7 in internal km / km/s units.
- Every Lambert `DV1` and `DV2` cached result is zero.
- All 35 Lambert `FinalState` position vectors are exactly identical:
  `[-34635399.02240507, -24067385.66544786, 3839.796195071856]` m.
  All 35 final velocity vectors are also identical:
  `[-742.0126620751473, -543.725488292007, -1630.695473128314]` m/s.
- F01 ends at UTC 2035-01-01 15:41:44.046993 with the above position. Coast2
  starts at that same UTC with position
  `[2475661.0267120823, 5493429.6880964153, 3519619.1504390808]` m.
  Their radii differ by tens of thousands of km; a coordinate-axis rotation
  cannot explain this discrepancy.
- Start and Coast1 also store different initial radii (609 versus 600 km
  altitude). These cached discrepancies do not prove the strategy cannot run:
  ATK may overwrite cached segment states during sequence execution.

The existence of start/end state fields alone therefore does not supply the
executed impulses. Even valid position endpoints and flight duration require
selecting a transfer branch. Impulse cost requires velocity immediately before
and after the SAME maneuver, expressed in the same frame; end velocity minus
start velocity across a coast includes gravity and is not impulse delta-V.

## Local reconstruction, not reproduction of the reported 9.9 result

`scripts/analyze_v3_tail27.m` reads XML with a DOM/XPath parser and fixes the
specified visit order, all waits, all flight durations, and the initial state.
It uses the existing zero-revolution Lambert/J2 transfer routine with no plane
penalty in branch ranking. Actual spacecraft arrival states are carried forward;
positions and velocities are never reset to targets or XML cached states.
The completed NEW pulse sequence then undergoes fixed-impulse independent J2
replay with all targets independently integrated.

Configured Start reconstruction:

- Raw delta-V: 15.277462865 km/s; 35 impulses.
- Independent nominal J2 visits: 35/35.
- Maximum target distance: 0.027885 km.
- Conservative altitude lower bound: 203.831962 km.
- Diagnostic elapsed time, including baseline state diagnostics: 9.111 s
  (excludes MATLAB startup and initial parsing/loading).
- ATK alignment remains unverified. This is neither the teammate's verified
  original trajectory nor a new cold-start search achievement.

A single sensitivity check using the cached Coast1 initial state (600 km), with
the same strategy, produced 15.285585264 km/s. The 0.008122399 km/s difference
rules out that particular 9 km initial-radius discrepancy as an explanation of
the roughly 5.3 km/s disagreement with the reported result in this local model.
It does not identify the actual cause. Executed pulse output, branch semantics,
dynamics/frame agreement, and the matching file/result version remain to check.

Artifacts, kept outside search directories:

- `runs/v3/diagnostics/tail27_20260925/comparison.mat`
- `runs/v3/diagnostics/tail27_20260925/teammate_plan.csv`
- `runs/v3/diagnostics/tail27_20260925/teammate_reconstructed_burns.csv`
- `runs/v3/diagnostics/tail27_20260925/cost_by_coverage.png`
- `runs/v3/diagnostics/tail27_20260925_cached_coast/comparison.mat`

Plan CSV columns: target ID, preceding wait, departure time, flight duration,
arrival time (seconds). Burn CSV columns: index, target ID (zero for baseline
rows), impulse norm, pre-burn speed, radius, inclination change in degrees,
full plane rotation in degrees, pre-burn a/e, post-burn a/e.

## Comparison with archived cold-start trajectories

Baseline A: `v3_policy_cold_dev_300_seed888_09`.
Baseline B: `v3_beam1_cs30_cold_dev_300_seed888_01`.
Their archived independent validation is retained as evidence; this diagnosis
replayed supplied witnesses for state statistics without a new search.

| Metric | Teammate strategy LOCAL reconstruction | Baseline A | Baseline B |
| --- | ---: | ---: | ---: |
| Raw total delta-V, km/s | 15.277463 | 18.339985 | 24.353595 |
| Impulses | 35 | 44 | 35 |
| Duration, days | 9.592126 | 9.736464 | 10.000000 |
| Cost through visit 8, km/s | 8.309517 | 4.835413 | 4.222105 |
| Cost through visit 27, km/s | 12.524470 | 11.746040 | 12.196960 |
| Remaining cost after visit 27, km/s | 2.752993 | 6.593945 | 12.156635 |
| Burns with inclination change >5 degrees | 0 | 0 | 3 |
| Maximum full plane rotation, degrees | 4.560 | 4.053 | 37.362 |
| Median impulse, km/s | 0.285339 | 0.382441 | 0.388461 |

Counts of burns over 60 s after the most recent supplied visit witness are
19, 6, and 1 respectively. This is a diagnostic of these schedules, not proof
that waiting causes lower cost or that all such baseline burns intentionally
implement delayed target transfers. Comparing equal visit counts also does not
equalize which targets remain or how difficult they are.

## Implications for the next decision

1. Stronger >5-degree inclination penalties cannot directly distinguish burns
   in baseline A: all its burns already fall below that threshold. Significant
   in-plane timing and orbital-energy costs remain. A nearly fixed inspector
   plane does not require matching each target plane: flybys can occur near
   intersections, without rendezvous velocity matching.
2. Both archived searches leave an expensive final set. Cheap early progress
   is not sufficient evidence of a cheap full mission. Remaining target
   geometry, time windows, and arrival velocity must influence construction.
3. Delayed targeted departure is worth prioritizing for investigation.
   `v3GuidedChildren` applies targeted burns at `node.t`; `v3Expand` samples
   waits for random burns but calls targeted time search at the original node.
   A delayed targeted burn can be expressed as a separate coast followed by a
   later expansion, but is not proposed as one combined guided action.
   At completion beam width one, the sole slot is taken before the zero-gain
   quota, leaving no guaranteed slot for the prerequisite coast.
4. This supports investigating a joint departure-wait/flight-time proposal
   and preserving useful coast opportunities within Beam P-ACO. It does not
   justify hard-coding the teammate order, fixed 45-degree inclination, or
   claiming that this change will achieve 6.1 km/s.
5. A reported 6 km/s competitor result is a useful benchmark, not evidence of
   its strategy. The invalid 6.76 km/s target-edge lower-bound claim cannot
   establish impossibility or attainable margin. A 9.9 km/s executed trajectory
   would be much stronger diagnostic evidence once its pulse record is aligned.

No search implementation was changed. The user's existing modification to
`v3CostCandidates.m` and the supplied ATK file were preserved and not committed.
The first diagnostic attempt stopped during XML parsing because a nested
helper reused the parent loop variable; this was fixed before propagation.
No failed attempt was counted as a physical result. No full regression suite
was run. Pre-change checkpoint: `89f3a82`, tag
`v3-rollback-20260925-tail27-pre`.
