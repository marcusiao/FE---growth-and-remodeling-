# Experimental q3 / BFGS candidate

Separate copy of the verified-ordering whole-eye kinematic adapter. Original q4/full-Newton model and results preserved. Same mesh, constitutive law, theta evolution, full tensor E1, BCs, pressures, schedules and compact output. q3 applies to volume and surface forms (8 points per hex). Reference-normal pressure only.

Inverse BFGS uses an exact tangent MUMPS factorization as its base, retained within a load/remodeling increment. Maximum10 updates before tangent reformation; reset every increment, on inadequate curvature, or on failed updated-direction line search. Maximum15 reformations and30 iterations. Low-rank two-loop updates avoid dense inverse storage. This follows FEBio's BFGS reuse principle, not an exact FEBio port: existing residual tolerance and Armijo residual line search retained; no FEBio cmax estimator or energy/displacement criteria adopted. Geometry checks unchanged at36 samples per cell, including trial iterates. Base linear residual checks apply to retained tangent solves, not to the updated quasi-Newton direction.

Timing records separate residual assembly, tangent assembly, determinant checks and base LU solves. bfgs_action includes base_linear_solve, so those two timings must not be added. LU timing includes initial factorization and subsequent triangular solves. Reformation/update counts saved per step.

No simulation authorized/launched for this candidate yet; full-run hold remains. Algebra and check-only results do not establish mechanical convergence, speedup or quadrature accuracy. Commands when explicitly requested:

```bash
python3 model_whole_eye_hex.py --check-only
python3 model_whole_eye_hex.py --smoke
```

FEBio references: https://febiosoftware.github.io/febio-docs/features/features/core_newtonstrategy_bfgs/ and https://febiosoftware.github.io/febio-docs/user/chapter3/3.3-control-section/ .

Verification: check-only results/20260924_143133_587929 PASS. check_algebra.py/algebra_report.json PASS: inverse BFGS matches dense formula within1.06e-15; PETSc MatLUFactorNum reports1 factorization for3 RHS solves. An initial obsolete reuse_factorization parameter failed and was removed; the installed solver reuses its unchanged matrix automatically. Setup source snapshot predates that API correction. Mechanical BFGS cycle not run.


## Selected full run — 2026-09-24
User selected this second setting after completed three-way comparison. Full100/500/100 launched in results/20260924_162255_916541 with persistent log launches/20260924_162255_668153/simulation.log. Earlier no-run/hold statements are historical and superseded for this launch. Short-cycle comparison passed; full trajectory is not yet complete or validated. Do not launch a duplicate.


2026-09-25: full run results/20260924_162255_916541 completed all700 steps in16h2m14s; saved strict residual and positive-sampled-determinant checks pass. Completion is not physical or mesh-convergence validation.


## CMT porous transfer handoff
See [CMT_POROUS_HEX_HANDOFF.md](CMT_POROUS_HEX_HANDOFF.md) for the mesh, numerical and loading setup to reuse in a separate porous candidate. This kinematic branch and its results remain preserved.

