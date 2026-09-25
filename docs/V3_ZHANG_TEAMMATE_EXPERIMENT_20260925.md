# Authorized teammate-order DP experiment

## Authorization and scope

The user explicitly requested adapting Zhang's open source code and testing
whether the teammate's visit order can give lower total delta-V. This is an
authorized historical-order comparison, not a cold-start run. No experiment
result is inserted into the cold-search work pool or elite archive. The
original reference code, source PDF, ATK file and target inputs are preserved.
The MATLAB implementation follows the project's language requirement.

## Source and adaptation

- Paper: https://arxiv.org/abs/2508.02904, especially Algorithm 1, Eq. (34),
  Section IV.E.3 and Table 1.
- Code: https://github.com/zhong-zh15/Multi_Flyby_Dynamic_Programming,
  revision 22126db5374f6eb465965a804d1a69706aefe1c0, MIT.
- Detailed theory/source boundary: V3_ZHANG_DP_SOURCE_CROSSCHECK_20260925.md.
- Ported core: epoch-pair dynamic programming, predecessor reconstruction,
  adaptive time tubes. Extension: retain all enumerated Lambert branches as
  arc states rather than preselect a rendezvous-minimum branch.
- Retain the corrected teammate reconstruction's initial orbit, first burn
  epoch and target order. Use target J2 ephemerides for proposal endpoints;
  two-body Lambert arcs and finite-arc altitude screening are proposals only.
- Remove GTOC launch allowances, terminal rendezvous, relative-speed limits,
  and year/day transfer constraints. A 1 s minimum positive arc duration is a
  recorded numerical approximation. No total-cost acceptance cap is applied
  in this diagnostic; physical costs above 6.1 are reportable here, not elites.
- The port follows the author's restricted visit-bound impulse topology.
  Delayed post-visit burns in the teammate's actual strategy are not preserved.
  Thus fixed-original-times versus DP-times isolates timing/branch changes
  within this restricted family; comparison with 9.987178114 km/s also changes
  topology and cannot isolate timing alone.
- Every improving proposal is repaired through the existing J2 shooting
  implementation, retaining actual propagated states. It is then checked by
  independent fixed-pulse replay of the spacecraft and all 35 targets. Repair
  may change the proposal branch; repaired costs are reported separately.

## Planned decision-bearing test

Nominal 480 s development cap, seed 888 (DP itself is deterministic), no ATK.
One fixed-time pass, then 9-point grids with half-width 21600 s and successively
halved widths around the previous best proposal. Source hashes, input digest,
configuration, all proposal/repair/verification outcomes and elapsed times are
stored in a new runs/v3/diagnostics/zhang_dp_tail27_* directory.

Question: Does coupled epoch/branch optimization in this restricted topology
produce an independently verified trajectory below the teammate baseline?
A failed or expensive proposal does not prove the order is intrinsically bad,
that free-maneuver DP is ineffective, or that additional budget would suffice.
No unrelated full regression suite or repeated optimization runs are planned.

## User extension and first execution

The user additionally requested an explicit delayed-maneuver variant. Added
`zhangdp.delayed`: layered search retains actual arrival velocity, propagates
each coast, chooses a departure wait and arrival epoch, and enumerates Lambert
branches. This is the general state-based stage recurrence with finite beam
pruning, drawing on the already reviewed Beam P-ACO limited-width construction
(https://arxiv.org/abs/1704.00702), not an exact two-epoch DP. It is not evidence
of implementing the entire Beam P-ACO algorithm in this comparison.

First executions used 300 s caps and ended earlier:

- visit_bound_01: 48.364 s, fixed-original-epochs 35/35 at
  30.772565429072 km/s; finest completed tube 35/35 at 13.194169707298 km/s.
- delayed_01: 14.159 s; fixed-original-epochs and waits reproduced 35/35 at
  9.987184175491 km/s, only 0.000006061380 km/s above the 9.987178114112
  reconstruction. This is reproduction, not an improvement.
- Several wider graphs failed because the public v3LambertBranches fallback
  forwarded an infinite revolution limit to a helper requiring an integer.
  Fixed the EXPERIMENT callers to compute the finite physical necessary bound
  for each endpoint/time pair first. No production search helper was changed.
- A three-label-per-epoch delayed beam lost its incumbent when adding options;
  nominal proposals could be worse. Before the necessary corrected comparison,
  expanded to 12 labels per epoch and explicitly protected the previous best
  proposal's complete chain. This extends the author's retained-center
  refinement to the beam case; it does not restore global optimality.
- Corrected comparison uses equal 300 s caps, same input order, initial orbit
  and first burn epoch. Both start with original epochs, then half-width
  21600 s and successive halvings. Visit-bound DP has 9 epochs per layer;
  delayed search has 5 epochs, up to 4 waits, 12 labels per epoch plus the
  protected incumbent. This is an engineering comparison, not isolated proof
  of the causal effect of delay under identical candidate evaluation counts.

Second executions: visit_bound_02 completed in 50.036 s with the same final
13.194169707298 km/s, all wider grids processed. delayed_02 completed in
44.398 s with 9.987180962379 km/s, not a meaningful improvement. Its nominal
costs still rose in some refinements, revealing an adaptation bug: the
optional-incumbent guard checked nargin < 11 although the function has ten
arguments, so the incumbent was always cleared. Corrected this to < 10 and
added a monotonic-incumbent assertion. A necessary delayed-only rerun will
check this specific correctness fix; the visit-bound run need not be repeated.

## Final measured comparison

All values below are RAW impulse-norm sums after full fixed-pulse independent
nominal-J2 verification, including separately integrated targets. They are
not penalty scores and not ATK execution results.

| Variant | Raw delta-V km/s | Visits | Maximum distance km | Altitude lower bound km | Duration days | Elapsed s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Teammate reconstructed baseline | 9.987178114112 | 35/35 | 0.027881900 | 203.831941426 | 9.592126224 | previously archived |
| Visit-bound DP, corrected _02 | 13.194169707298 | 35/35 | 0.039492263 | 241.257473607 | 9.588219974 | 50.0356583 |
| Delayed finite beam, corrected _03 | 9.986974010255 | 35/35 | 0.037421303 | 203.832487448 | 9.592126224 | 55.2964191 |

The delayed variant improves the reconstructed baseline by
0.000204103856 km/s = 0.204103856 m/s = 0.002043659%. This is a small local
improvement, not a competitive-cost breakthrough. The best delayed result
comes from iteration 5; its two-body proposal was 9.980409149728 km/s, which
MUST NOT replace the verified 9.986974010255 km/s in reporting.

The final delayed candidate keeps all original planned arrival epochs. Its
substantial timing changes postpone burns 29 and 30 by 675 s each. Other
departure-time differences are below 0.04 s and arise from the reconstructed
plan/wait values. Downstream impulses change as actual velocities propagate,
so local impulse changes must not be mistaken for net mission improvement.
It retains 19 post-visit waits longer than 60 s and 35 impulses.

The corrected delayed implementation retains its incumbent through every
layer and checks that the final graph cost cannot increase beyond 1e-8 km/s.
All eight iterations of _03 completed without the preceding errors. The
source hashes were unchanged throughout each final run; the summary utility
also checks recorded independent-pass status, the raw pulse sum, ordered
witnesses and nonnegative waits from the saved schedules.

Both final runs used 300 s soft caps but finished their eight configured
resolutions early. These are fixed-iteration development experiments, not
300 s exhaustive searches or convergence proofs. Including first/debug runs,
measured algorithm time was about 212 s, excluding MATLAB startups. No formal
24-minute batch, ATK session or unrelated regression suite was run.

## Interpretation

- Visit-bound time/branch optimization reduces its own fixed-time baseline
  from 30.772565429072 to 13.194169707298 km/s. It does not beat the original
  delayed strategy.
- Preserving delays reproduces the 9.987 km/s strategy, and optimizing wait
  choices yields the small independently verified improvement above.
- This supports keeping delayed maneuvers in the construction family. It
  does NOT quantify their isolated optimal benefit: the two algorithms have
  different state spaces, candidate grids and pruning rules.
- The delayed beam still prioritizes prefix cost and has finite width. A
  locally more expensive arrival velocity with a better future may be lost.
  Keeping one incumbent prevents regression but does not fix this exploration
  limitation. This remains a candidate explanation, not a measured cause.
- The initial orbit, first burn epoch, order and time-tube centers use the
  user-authorized teammate reference. No cold-start claim or acceptance-line
  success follows. No search-model change or production integration occurred.

## Artifacts and reproduction

- `scripts/test_zhang_dp_tail27.m`: visit-bound MATLAB adaptation.
- `scripts/test_zhang_delayed_tail27.m`: delayed finite-width extension.
- `scripts/+zhangdp/`: standalone recurrence implementations and MIT notice.
- `runs/v3/diagnostics/zhang_dp_tail27_visit_bound_02/experiment.mat`.
- `runs/v3/diagnostics/zhang_delayed_tail27_delayed_03/experiment.mat` and
  `best_new.mat`; the latter contains the new schedule and independent report.
- `runs/v3/diagnostics/zhang_tail27_comparison_03/`: summary CSV/MAT,
  three pulse-table CSVs and comparison.png, visually inspected.
- Both experiment MAT files retain the seed, configuration, model/source
  hashes, input digest, intermediate proposals, failures and actual times.
- The pre-existing uncommitted change in v3CostCandidates.m was left intact;
  it is included in the implementation signature but not part of this port.

From the simulation directory, MATLAB commands for new output directories:

```matlab
addpath('scripts');
visitBound = test_zhang_dp_tail27(300);
delayed = test_zhang_delayed_tail27(300);
delayed.best_new.verification
delayed.improvement_km_s
```

These commands intentionally load the teammate diagnostic as an explicitly
authorized comparison input. They must not be substituted for cold-start
experiment entry points. Existing output directories are never overwritten.

Status: both requested variants implemented and compared; all experiment
processes ended. No further search batch is running.
