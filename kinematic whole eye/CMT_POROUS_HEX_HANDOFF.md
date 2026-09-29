# Whole-eye hex handoff: kinematic baseline to CMT porous constitutive model

Prepared 2026-09-25 from the current scripts, mesh reports and completed run. This is a handoff, not an implemented CMT porous model. User wants the same whole-eye mesh, loads and constraints with theta acting through the CMT porous strain-energy density (SED, also called SDE in discussion), rather than through Fg. Preserve this kinematic branch and its results; develop the porous adaptation in a separate candidate folder. Read project-root README_CODEX.md and PROJECT_FILE_INDEX.md first. No new simulation was launched for this document.

## 1. Authoritative baseline and environment

Project root (Windows):
`C:/Users/E1521285/Desktop/Linux_shared_file/ONH_FEniCS_Project`

Project root (WSL):
`/mnt/c/Users/E1521285/Desktop/Linux_shared_file/ONH_FEniCS_Project`

All paths below are project-relative unless otherwise stated.

- Selected numerical candidate: `05_main_model/kinematic_whole_eye_hex/q3_bfgs/`.
- `model_whole_eye_hex.py`: mesh/tag import, material assignment, weak forms, full tensor stimulus, evolution, cycle and output.
- `hex_geometry.py`: source-order hex interpolation, integration and determinant sampling.
- `bfgs_solver.py`: selected strict BFGS/MUMPS implementation.
- Completed full run: `05_main_model/kinematic_whole_eye_hex/q3_bfgs/results/20260924_162255_916541/`.
- Executed source snapshot: that run's `source/`; hashes and mesh manifests in `summary.json`. Prefer this snapshot when reproducing the exact baseline.
- Persistent full-run log: `05_main_model/kinematic_whole_eye_hex/q3_bfgs/launches/20260924_162255_668153/simulation.log`.
- Existing Ubuntu-24.04 WSL, legacy DOLFIN2019.2.0.64.dev0, serial, MUMPS. No FEniCSx migration or software installation is needed for this mesh adapter.

The selected setting is q3/BFGS with original strict residual criteria and36-point checks. The sibling `q3_febio_controls/` experiment is NOT the user's selected baseline.

## 2. Mesh input and essential legacy-DOLFIN ordering fix

Original, read-only source bundle:
`01_mesh_conversion/whole_eye_hex_no_bm_20260923/mesh/`

Actual runtime bundle to import:
`05_main_model/kinematic_whole_eye_hex/reordered_mesh/`

Both bundles contain `whole_eye_domain.xdmf/.h5` and `whole_eye_facet_region.xdmf/.h5`. Keep each pair together. Original `mesh_metadata.npz` and `conversion_report.json` supply mapping/identity checks; runtime `verification.json` supplies runtime-file hashes. Use the reordered bundle for DOLFIN, not the original XDMF directly.

Mesh: **36,864 hex8 cells,42,494 vertices**. BM shells removed; shared solid interfaces remain connected. No BM stiffness, fluid mechanics or fibre mechanics is implemented. Source fibre/frame metadata exists but is not used by this isotropic baseline.

Original import failed with `Cell is not orderable`. Applying only a local VTK-to-UFC permutation did not solve the global-order requirement. `analyze_ordering.py` in the parent directory constructed globally consistent vertex numbering using149 parallel-edge classes and integer vertex potentials. Every local reordered hex was verified against the48 cube symmetries. `audit_ordering.py` checked node bijection and preservation of each cell's edges/faces. Coordinates were only permuted, never moved; no remeshing or tetrahedral subdivision occurred. `verify_ordered_mesh.py` exported and round-tripped the separate runtime bundle.

Read `05_main_model/kinematic_whole_eye_hex/ORDERING_FIX_20260924.md` and runtime `ordering_audit.json`/`verification.json`. This construction succeeded on this mesh; it is not a generic guarantee for all hex meshes. Do not simply sort arbitrary eight-node connectivity or regenerate the already verified bundle unnecessarily.

`Model.__init__` matches runtime vertices by coordinates and cells by vertex membership to the original metadata, then verifies material and boundary tags. Sorted node tuples are matching keys, not interpolation connectivity. Preserve `source_order_cells`: NumPy geometry calculations require VTK source order, while DOLFIN uses its own ordered connectivity. Never assume cell/vector array order matches source order.

## 3. Tissue IDs and baseline material values

Read tissue tags from XDMF attribute **`name_to_read`**, dimension3. These are material IDs, not FEB block `domain_id`.

| Tissue tag | Tissue | Hex count | E (MPa) | Status |
|---|---|---:|---:|---|
| 1 | Sclera | 13,536 | 3.0 | Reused from previous kinematic ONH model |
| 2 | Choroid | 4,824 | 0.2 | Reused |
| 4 | Retina | 6,696 | 0.1 | User-authorized exploratory placeholder |
| 5 | Optic nerve (ON) | 4,896 | 0.1 | Placeholder |
| 6 | Pia | 504 | 0.3 | Placeholder |
| 7 | Dura | 1,080 | 3.0 | Placeholder |
| 8 | Lamina cribrosa (LC) | 1,836 | 0.3 | Reused |
| 9 | Prelaminar tissue (PT) | 3,492 | 0.1 | Reused; ONLY remodeling tissue |

All baseline Poisson ratios are0.49. `mu=E/[2(1+nu)]`, `lambda=E*nu/[(1+nu)(1-2nu)]`. Coordinates use mm, stress/pressure MPa, force N, volume mm3. E and nu are baseline parameters; converting them into porous intrinsic/effective moduli must follow the selected porous law, not an invented equivalence.

Tissue tag3 is absent. **Boundary tag3 remains valid IOP.** In the older tetrahedral FEB kinematic model PT was tag3; here every PT mask must use9. Do not carry old tissue IDs into the new model.

Reference PT volume:0.92450463615135mm3; LC:0.38400851071697mm3. Full per-tissue volumes are in the completed run's `summary.json` setup block.

## 4. Boundary IDs and exact applied conditions

Read facet tags from `name_to_read`, dimension2. Markers do not impose BCs automatically.

| Facet tag | Source name | Facet count | Baseline condition |
|---|---|---:|---|
| 1 | ZeroDisplacement2 | 648 | ux=uy=uz=0 |
| 2 | symmetry | 2,128 | **uz=0 ONLY**; ux/uy free unless also on a clamped boundary |
| 3 | IOP | 2,772 | Inward reference-normal traction, maximum0.0047MPa |
| 4 | CSFP | 756 | Inward reference-normal traction, maximum0.00172MPa |
| 5 | ZeroDisplacement3 | 756 | ux=uy=uz=0 |

The assembled constraints contain6,960 unique displacement DOFs. Unloaded exterior faces have zero natural traction; shared interior faces are not exterior loads. No imposed scleral canal radial displacement, no LC spring/foundation and no additional bottom support are present.

**Pressure convention matters:** outward undeformed normal N, traction t=-p*N per reference area. With residual internal minus external, pressure terms are:

```python
R += factor * (0.0047*dot(N,v)*ds(3) + 0.00172*dot(N,v)*ds(4))
```

This is dead reference-normal traction, NOT follower pressure on the deformed surface. The selected CLI permits `--pressure-mode reference` only. To reproduce the same simulation conditions retain this convention, even if the user's FEBio run uses another convention.

## 5. What was changed from tetrahedra to hexes

- All displacement elements are continuous trilinear **Q1**, implemented as `VectorFunctionSpace(mesh,'Q',1)`;127,482 displacement DOFs. There is no partial quadratic region and no displacement-pressure mixed formulation.
- Cell scalar fields use `FunctionSpace(mesh,'DQ',0)`; full strain fields use `TensorFunctionSpace(mesh,'DQ',0,shape=(3,3))`. Explicit `qmap`, `tmap` and vertex-to-DOF maps are retained. Do not index raw function vectors by cell ID without these maps.
- Hex gradients vary within cells. The old linear-tetrahedron constant-gradient, one-J-per-cell and edge-determinant/6 shortcuts were replaced by trilinear hex geometry and FE integration.
- `F=I+grad(u)`, `E=0.5*(F.T*F-I)` includes all9 tensor components. A local FE projection obtains the cell-volume-averaged tensor, symmetrized and diagonalized with `eigvalsh`; E1 is its largest eigenvalue. It is NOT the average of pointwise largest eigenvalues, and not an axial or diagonal-only approximation.
- FE volume/surface measures use **quadrature degree3**:8 volume integration points per hex in this installed environment. Degree2 also gives8; degree4 gives27.
- Independently, `hex_geometry.volumes` uses3x3x3 geometry integration. `sampled_J` checks27 Gauss locations +8 corners +centre =36 samples, computing current/reference geometric Jacobian determinant ratio. This is separate from q3 force integration and is not a whole-cell positivity proof.
- PT/LC geometric volume uses `assemble(J*dx(tag))`; visualization J is a cell average. Retain integrated volume logs rather than relying on ParaView appearance.

## 6. Constitutive replacement boundary: kinematic versus CMT porous

Current kinematic law:

```text
Fg = theta^(1/3) I
Fe = F inv(Fg)
Je = det(Fe) = J/theta
W = mu/2*(tr(Fe.T Fe)-3) - mu*ln(Je) + lambda/2*ln(Je)^2
```

Energy is per original reference volume, with no theta prefactor or Jg multiplier. Here theta is a **natural-volume ratio**. This is why theta<1 changes the unloaded preferred geometry.

Intended CMT porous replacement: use the explicitly selected existing porous `W(F,theta,...)` and its porosity/solid-compressibility equations. Theta should represent remaining structural content according to that law. Do not retain Fg shrinkage in addition to porous degradation. Do not assume porosity equals1-theta unless the selected law, initial porosity and volume convention establish that relation. The exact porous version/parameters are not selected or reproduced by this handoff.

Integration points in the code to adapt:

1. Replace the Fg/Fe/Je and W construction in `Model.__init__`; keep F, total-strain E, mesh and loading infrastructure.
2. Derive `P=d.diff(W,F)` and the residual/tangent from the actual new energy, with theta frozen inside each equilibrium solve. Verify that W and its derivatives use the intended reference-volume convention.
3. Remove/rename misleading `Theta_Jg` output; replace `sampled_Je_min=J/theta` with diagnostics appropriate to the porous law. Without Fg that quotient is not automatically an elastic volume ratio.
4. Add the selected law's porosity, pore-closure and intrinsic-solid diagnostics and admissibility checks. Positive J alone may be insufficient; check logarithm/denominator domains using the actual law. Do not reintroduce an incompressible-solid J-s barrier if the selected law allows intrinsic solid compression after closure.
5. Retain full-tensor E1 extraction, but audit whether the intended porous stimulus/update is identical to this baseline. Do not silently copy the same symbol theta while changing its meaning or add theta scaling twice.

Project history points to porous experiments under `05_main_model/validation/porous_onh_pt_lc/`, including `lc_loaded_foundation/`; read the exact selected script, README and verification reports before extracting its law. These are read-only references. Some historical porous formulations have documented material-instability/strong-ellipticity concerns. Kinematic convergence does not validate them. In particular, **do not import the LC foundation, old tetrahedral material/facet IDs or partial-CG2 machinery** into a same-BC whole-eye comparison.

## 7. Baseline evolution and phase schedule

Initial theta=1 everywhere. Only PT9 updates. Baseline E1crit=.005:

```text
overload = clip((E1-E1crit)/E1crit, 0, 5)
theta_next = clip(theta - .01*((theta-.5)/.5)^2*overload*.1, .5, 1)
```

E1 comes from the previous converged state; theta is updated once before each remodeling equilibrium solve. Non-PT theta remains1. This describes the executed kinematic baseline, not an instruction to override the selected porous model's evolution without review.

- Loading100 steps: both pressures ramp together by factor step/100; theta fixed.
- Remodeling500 steps: both pressures held at full magnitude; PT evolves then equilibrates each step.
- Unloading100 steps: both pressures reduce together by factor1-step/100; theta frozen at end-remodeling state.
- Smoke schedule10/3/10, same endpoint pressures. Remodeling dt=.1 is an evolution parameter, not a calibrated clinical time unit.

## 8. Selected solver and output workflow

Use strict `q3_bfgs/bfgs_solver.py` as the numerical reference: exact tangent/MUMPS factorization at each increment, inverse BFGS updates, max10updates before reformation, curvature/line-search safeguards, max15reformations/30iterations. Reuse within increments only; no cross-increment factor reuse. Stop when residual L2 <=max(1e-9,1e-8*initial_residual). Up to24halved line-search trials, Armijo coefficient1e-4; positive finite sampled J>1e-8 before trial residuals. Base LU relative residual <=1e-6. These settings differ from FEBio's displacement/energy stopping criteria.

BFGS convergence and stability must be retested for the selected porous energy; it may have different curvature. Do not suppress convergence failures or alter constitutive parameters automatically. Timing `bfgs_action` includes `base_linear_solve`; do not add them twice.

Each phase writes separate XDMF/H5 with its start, first step, every10th step and endpoint. Every step records scalar diagnostics in `measurements.csv` and `summary.json`. Save source snapshots, hashes and input manifests. Existing outputs are not restart checkpoints. Use unique result folders; preserve failed runs.

When relocating the adapter, deliberately repair `ROOT`, `MESH` and `RUNTIME_MESH` paths; their current parent-depth expressions assume the q3_bfgs location. Reference the existing verified runtime mesh read-only. The kinematic `launch_full.py` checks hashes against its own short run and is not a generic launcher for a new constitutive model.

## 9. Evidence and next verification gates

Completed kinematic full cycle:100/500/100,57734.290s (16h2m14s), status `COMPLETED_EXPLORATORY_NOT_VALIDATED`. All saved residuals pass, sampled minima J=.658754817 and Je=.884797793. User reports expected qualitative behavior. This is not mesh/locking/biological validation.

Supporting reports: parent `ORDERING_FIX_20260924.md`; parent `verification/ordering_integration_20260924_091406_813610/report.json` (original q4 prescribed-state tangent/normal/output tests); q3 `results/20260924_143133_587929/summary.json` setup check; `algebra_report.json` BFGS/LU reuse test; q3 `results/20260924_143754_479346/COMPARISON.md` short comparison. Original q4/Newton short run took6578.875s, strict q3/BFGS2023.373s. Numerical gains are not a porous runtime guarantee.

For the new porous candidate, verify mesh/tag counts,6,960BC DOFs, reference volumes, pressure signs and full-tensor affine extraction first. Verify zero-load reference stress at the selected initial porosity/content, finite-difference energy/residual/tangent consistency at heterogeneous theta, closure transition and constitutive domains. Then run an explicitly authorized short cycle and compare loading, PT/LC volumes, displacement and content evolution before a full cycle. A changed porous healthy-state law may produce a different loading response despite identical E/nu labels; document it. Preserve q3 and strict checks as the starting numerical baseline unless a separately recorded test justifies changes.

Unresolved: displacement-only locking at nu=.49, spatial/quadrature convergence, placeholders for four tissues, stability/calibration of the porous law, whole-cell deformed positivity and clinical validity. Keep numerical adaptation separate from constitutive choices.
