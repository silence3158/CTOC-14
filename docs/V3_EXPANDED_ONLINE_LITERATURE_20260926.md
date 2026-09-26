# Expanded online research: trajectory-centered multi-target flybys

Date: 2026-09-26. Scope: literature research requested by the user, not an
implementation round or an optimization experiment. No solver changes or runs.

## Question

Who has designed trajectories that naturally pass multiple targets, rather than
constructing every maneuver as an independently targeted transfer? How can the
user's proposal be expressed more precisely?

New online searches used Crossref, publisher pages, Research Square, J-STAGE,
OpenAlex, and reference tracing. Previously discussed Zhu, Zhang (2025), and
Huang (2025) papers are not counted as new findings. Evidence levels below are
deliberately different: metadata, abstract, and method-section checks.

## 1. Trajectory propagation discovers candidate visits

Quan Jing, Zhixin Hao, Mingtao Li (2024), "Trajectory optimization of flybys of
multiple irregular satellites of Jupiter with Galilean moons gravity assist."

- Published DOI: https://doi.org/10.1007/s10509-024-04305-7
- Public manuscript: https://www.researchsquare.com/article/rs-3892642/v1
- Evidence: downloaded public PDF; checked introduction, sections 3.1-3.3 and
  4.4, and references. Published-version link confirmed on Research Square.
- Section 3.2 propagates the current orbit without a new maneuver, compares
  spacecraft and target positions at matching epochs, and retains nearby
  target/time opportunities. Lambert transfers then correct these opportunities;
  beam search retains branches. It includes cases without gravity assists.
- Actual relevance: trajectory-led candidate discovery and small corrections.
  This is NOT a demonstration of jointly fitting many targets to one coast arc:
  subsequent construction still selects individual targets and Lambert legs.
- Section 4.4 retains proportions of different branch classes because selecting
  only by immediate delta-v can eliminate useful alternatives.
- Jupiter-specific numerical costs and gravity-assist benefits do not establish
  CTOC14 performance. No public implementation was verified in this review.

## 2. One impulse, two or three shared encounters

Cunyan Xia, Gang Zhang, Yunhai Geng (2022), "Coplanar multi-target interception
with a single impulse" (see publisher for the original Chinese title).

- https://hkxb.buaa.edu.cn/CN/10.7527/S1000-6893.2021.25093
- Journal: Acta Aeronautica et Astronautica Sinica, 43(3), 325093.
- Evidence: publisher Chinese and English abstracts and reference list checked.
- Two-target case: Gibbs three-position orbit determination reduces the problem
  to nonlinear equations in two free variables when a maneuver or encounter epoch
  is specified. Newton-Raphson solves these equations; optimizing the maneuver
  epoch obtains a fuel-optimal solution within the studied setting.
- Three-target case: Lambert-based reduction to two free variables; numerical
  root solving and Pork-Chop initialization.
- Actual relevance: several encounters jointly determine ONE ballistic transfer
  orbit. This directly supports a multi-target arc construction primitive.
- Studied coplanar setting; no claim here of a general 35-target J2 solution or
  of globally optimal arbitrary multi-impulse schedules.

Haoxiang Su, Zhenghong Dong, Lihao Liu, Lurui Xia (2022), "Numerical Solution for
the Single-Impulse Flyby Co-Orbital Spacecraft Problem."

- https://doi.org/10.3390/aerospace9070374
- Public PDF: https://mdpi-res.com/d_attachment/aerospace/aerospace-09-00374/article_deploy/aerospace-09-00374.pdf
- Evidence: downloaded PDF; checked sections 2-3, conclusions, and references.
- One maneuver takes a chaser onto an orbit that encounters two targets sharing
  the same orbit at different phases. Six unknowns are maneuver epoch, three
  velocity components, and two encounter epochs. Geometric relations reduce the
  problem to one variable at a specified maneuver epoch.
- Actual relevance: the first encounter does not trigger another maneuver;
  departure design simultaneously serves both encounters.
- Its co-orbital two-body reduction is specialized. Reported numerical speedups
  and example costs were not independently reproduced here.

## 3. Ballistic building blocks and global sequence design

Giuseppe Cataldi, Salvo Marcuccio (2022), "Design of 3-D trajectory sequences for
multiple asteroid flyby missions."

- https://doi.org/10.1007/s42401-022-00166-6
- Evidence: Crossref publisher-deposited abstract and authors checked; publisher
  full text could not be retrieved in this session.
- Abstract explicitly describes impulsive maneuvers connecting ballistic coast
  arcs and a deterministic building-block approach for multiple Apollo asteroid
  flybys, considering target count and propellant consumption.
- Relevance: coast-arc composition as a trajectory-design representation.
- Abstract alone does NOT establish that each building block visits multiple
  targets, or that all connections are optimized jointly. Do not claim either.
- Earlier related article: "A 2-D Trajectory Design Algorithm for Multiple
  Asteroid Flyby Missions" (2020), https://doi.org/10.1007/s42496-020-00067-x .
  Metadata checked only; repository landing page was inaccessible.

Mai Bando, Hiroshi Yamakawa (2010), "Orbital Design for Multiple Flyby Missions."

- https://doi.org/10.2322/tastj.8.Pd_9
- Publisher: https://www.jstage.jst.go.jp/article/tastj/8/ists27/8_ists27_Pd_9/_article
- Evidence: publisher abstract checked.
- Uses canonical-system generating functions to evaluate two-point boundary
  problems and select fuel-minimizing multiple-flyby sequences. Demonstration
  uses Hill dynamics and prescribed inter-encounter intervals.
- Relevance: a longer-established multi-encounter optimization reference, but
  less directly aligned with free-time natural-arc discovery than items 1-2.

## Further leads, not established methodological recommendations

- Xia, Zhang, Geng (2021), "Two-target interception problem with a single
  impulse": https://doi.org/10.1016/j.ast.2021.107110 . Title/authors verified
  through Crossref; Su's paper describes its noncoplanar treatment. Original
  method not independently checked here.
- Bellome et al. (2024), "Modified dynamic programming for asteroids belt
  exploration": https://doi.org/10.1016/j.actaastro.2023.11.018 . Metadata and
  citation in Jing verified; original abstract/full text not obtained.
- Russell et al. (2023), "Global trajectory optimization, pathfinding, and
  scheduling for a multi-flyby, multi-spacecraft mission":
  https://doi.org/10.1016/j.actaastro.2022.07.007 . Metadata only.

## Synthesis for discussion, not a published algorithm claim

The user's idea is most precisely a joint design of a controlled trajectory and
its encounter events. Let q contain the initial state and the complete impulse
schedule, and r(t;q) be the resulting propagated trajectory. Each target's
encounter time is free and its visit is expressed by

    min over t in [0,T] ||r(t;q) - r_i(t)|| <= epsilon.

Minimize the sum of impulse-vector norms. Encounter order follows the selected
event times. This formulation does not require a distinct impulse per target.
For numerical optimization, free witness times can replace the inner minimum;
multiple approach windows remain a nonconvex search issue.

The three useful concepts are (a) propagate first to discover opportunities,
(b) fit several encounters jointly to a shared natural arc, and (c) evaluate and
adjust connections using the entire physical trajectory. Items 1 and 2 provide
particularly concrete new literature for the first two concepts. This review
does not establish one paper implementing all three for CTOC14.

No assertion is made that these methods attain 6.1 or 8 km/s on this instance.
No new local cost or feasibility evidence was generated. Some services returned
429 or access challenges; these limitations are reflected in evidence levels.
