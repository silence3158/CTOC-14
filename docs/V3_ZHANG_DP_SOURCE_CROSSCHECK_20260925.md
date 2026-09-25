# Zhang 2025: paper and source crosscheck

Date: 2026-09-25. Scope: read-only research requested after the user paused
algorithm development. No search, reference executable, MATLAB benchmark, or
ATK run was performed in this review. Search sources and original inputs were
not edited. This note is separate from V3_ZHANG_DP_REVIEW_20260925.md, which
appeared in the shared worktree during this review and was left untouched.

## Sources and identity

- Zhong Zhang, Xiang Guo, Di Wu, Hexi Baoyin, Junfeng Li, Francesco Topputo,
  *Global Optimality in Multi-Flyby Asteroid Trajectory Optimization: Theory
  and Application Techniques*, arXiv:2508.02904v1, 2025-08-04.
  https://arxiv.org/abs/2508.02904
- Local PDF: `literature/Zhang et al. (2025) .pdf`.
  SHA-256: `2FC12BF13B824E428E5CDE3D5E68054351F7DA99988599FEF3FDED4DA9722EEE`.
- Source: https://github.com/zhong-zh15/Multi_Flyby_Dynamic_Programming
- Local repository HEAD and GitHub HEAD checked:
  `22126db5374f6eb465965a804d1a69706aefe1c0`.
- DP_flyby_2impulse.cpp, GTOC4_Problem.h, and Lambert.cpp Git blob hashes
  match that upstream revision. The local repository has modifications to
  DP_flyby_2impulse_user_function.cpp and main.cpp; do not describe the whole
  checkout as byte-identical to upstream. Findings below refer to local code.
- Read paper abstract, Sections II-IV, Algorithm 1, Table 1, case descriptions
  and conclusion; also rendered and inspected Figure 4 on PDF page 19.
- The README calls the paper submitted in 2024. The supplied document is the
  2025 arXiv version; journal publication status was not established here.

## Main idea

The paper solves trajectory optimization UNDER A GIVEN VISIT SEQUENCE.
It does not prove that low-cost tours are stationkeeping trajectories. Its
relevant contribution is joint time/arrival-state optimization using Bellman's
principle, with a specialized low-dimensional ballistic-arc implementation.

At a flyby, spacecraft arrival velocity is generally different from the
target velocity. With impulses at visit epochs, the connection cost is

    dv_k = norm(v_departure_next_arc - v_arrival_previous_arc).

Neither endpoint is reset to target velocity. A currently expensive incoming
arc can make the next connection cheaper; selecting the cheapest individual
leg is therefore insufficient.

For fixed target order, endpoints exactly at target positions, and a chosen
ballistic branch, two consecutive epochs determine the incoming velocity.
The special-case state can therefore be `(t_previous, t_current, branch)`.
The supplied code deterministically preselects one branch per epoch pair and
omits the branch index. It retains the cheapest history for each resulting
epoch pair after evaluating the incoming/outgoing velocity difference.

A branch-aware recurrence for this RESTRICTED construction family is:

    F[k+1](t[k], t[k+1], b_new) = min over t[k-1], b_old of
        F[k](t[k-1], t[k], b_old)
        + norm(v_out(b_new) - v_in(b_old)).

The launch contribution and terminal contribution must follow the actual
problem. CTOC14 has no mandatory terminal velocity-match impulse.

General free-maneuver trajectories cannot be compressed to two epochs alone.
The full Markov state must retain actual r/v, absolute time, visited targets,
and any other state needed for future constraints. At a <=1 km encounter,
the position offset is also not uniquely specified by target and epoch.

## Code evidence

Paths below are relative to `reference code/Multi_Flyby_Dynamic_Programming`.

1. `code/src/DP_flyby_2impulse.cpp:60-107`: for each next epoch, group
   predecessors by current epoch; evaluate transitions and keep the cheapest
   resulting state for that ordered epoch pair. The retained rv is the
   selected arc's arrival state. This is not one winner per target alone.
2. `code/src/DP_flyby_2impulse_user_function.cpp`, compute_rv_dv:
   subtract last_state.rv velocity from new departure velocity. The GTOC4
   path adds a terminal rendezvous burn and subtracts a 4000 m/s launch
   allowance; neither belongs in CTOC14. GTOC11 also uses a 2000 m/s relative
   flyby-speed constraint and a 6000 m/s launch allowance.
3. `code/src/GTOC4_Problem.h:60-75` calls solve_multi using TARGET states
   at both endpoints. `code/include/Lambert.cpp:338-365` enumerates revolution
   counts and left/right solutions, but keeps only the smallest sum of two
   rendezvous-like endpoint impulse magnitudes. DP then uses spacecraft
   velocity continuity, a different criterion. It is exact at most for the
   resulting restricted graph, not for all physical Lambert branches.
4. `code/include/Lambert.cpp:343-355`: the check of Lambert flag0 is commented
   out; candidate velocities are used in cost evaluation and the wrapper can
   set its own success flag. This is a static numerical-validation risk,
   not evidence that the paper's published trajectories failed. A port needs
   solver-status, finite-state, endpoint-residual, and physical checks.
5. `code/src/DP_flyby_2impulse.cpp:389-398`: recursive refinement halves
   time step AND window around the previous solution, using a cost threshold
   total_dv + 1e-3. The GTOC4 units are m/s, so this margin is 0.001 m/s,
   NOT 1 m/s. The previous epochs occur at the centers of the new grids.
6. `code/src/main.cpp:27-40` loads a supplied GTOC4 sequence. The initial
   epoch windows use the midpoint of its overall time span and stage-wise
   uniformly spaced centers. The line using EACH original visit epoch as
   the center is commented out in possible_t_values. This is still a
   reference-derived setup, but not per-visit centering on the champion.
7. check_t restricts GTOC4 inter-visit flight times to 1 day through 1 year.
   These are not CTOC14 restrictions. The implementation uses heliocentric
   two-body arcs, not the project's Earth-centered J2 dynamics.

## What the paper demonstrates, and what it does not

Table 1 reports that a 32-day fixed time grid gave 228,369.5 m/s in 0.311 s,
while fixed plus adaptive refinement gave 25,210.4 m/s in 12.428 s for its
GTOC4 impulsive example. This is evidence in the AUTHOR'S problem that time
resolution and coupled refinement can matter dramatically. It is not our
runtime measurement and cannot predict CTOC14 improvement or runtime.

Section IV.E.3 explicitly says that shrinking the search to a tube around
one previous solution sacrifices the theoretical global-optimality guarantee.
An exact DP optimum on a finite graph is distinct from a global optimum of
the continuous free-maneuver problem, and from an optimum over target orders.

The paper's Eq. (33) relates error to N times a maximum stage error under its
formulation. We have not established the relevant error bounds, feasible-set
coverage, or continuity assumptions for CTOC14. A few two-body/J2 comparisons
or a dimensional acceleration estimate would not establish such a bound.
We therefore cannot use this formula as a numerical optimality certificate.

## Transfer to V3: proposed direction, not implemented

- Use coupled epoch AND branch comparison to construct trajectories from
  candidate sequences generated in the current cold-start run. Keep outer
  sequence/initial-orbit exploration and the free-maneuver optimizer.
- Treat visit-bound impulses as one restricted construction family, not a
  replacement for delayed burns, zero-visit transitions and multi-target arcs.
  The teammate strategy contains delayed post-visit departures, which would
  be excluded by a direct copy of the specialized two-epoch DP.
- Preserve multiple different coarse paths for refinement when practical;
  this mitigates dependence on one initial tube, but is our proposed
  extension and does not restore a global guarantee.
- In the restricted graph, compute a candidate arc once per epoch pair and
  branch, then reuse its endpoint velocities across predecessor-cost
  comparisons. The reference loops repeatedly call compute_rv_dv for the same
  pair. State matching, branch retention and J2 validity still govern reuse.
- Any two-body surrogate is a proposal generator. J2 shooting, full-path
  altitude checks and independent fixed-pulse 35-target verification remain
  necessary. No quantitative claim of small J2 correction cost is supported
  by this review, especially across low-perigee long arcs.

The raw cost of ANY independently verified feasible restricted-family
trajectory is an upper bound on the unrestricted optimum (minimization).
It is not a lower bound. An unverified two-body surrogate has neither status
for the J2 task. A high-cost restricted DP outcome cannot prove that another
order, branch family, time grid, or free-maneuver topology is also expensive.

The review does not resume algorithm development or authorize historical
solutions as cold-search inputs. No new result meets the user's acceptance
criterion as a consequence of this reading.
