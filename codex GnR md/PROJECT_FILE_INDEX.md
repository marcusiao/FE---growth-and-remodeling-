# Project file index
Read [README_CODEX.md](README_CODEX.md) first. Paths below are relative to this workspace. Every pre-existing research script/report/result is **reference/read-only by project policy**; no filesystem permission changes were applied. New development must use a separately named candidate in the user-specified folder. A copy is not a new scientific version.

## Main model and full-tensor work
| Path | Purpose and status | Supersession / evidence |
|---|---|---|
| [05_main_model/reference/E1_full_tensor_clean.py](05_main_model/reference/E1_full_tensor_clean.py) | Main current structural-weakening-only E1 research reference; full 3x3 eigvalsh, CG1 displacement. Read-only. Sibling domain/facet XDMF/H5 included. | Clean derivative of diagnostic implementation; identical hash to locking reference. Strain calculation preferred, CG1 mechanics still subject to locking. |
| [05_main_model/current_candidate/README.md](05_main_model/current_candidate/README.md) | Active PT+LC regional candidate, deferred buffered variant, and preserved coarse-CG2 history; production unvalidated. | User selection and q2 loading test supersede earlier buffered preference; historical regional sensitivity remains recorded. |
| [02_full_tensor/implementation/project/E1_full_tensor.py](02_full_tensor/implementation/project/E1_full_tensor.py) | Full-tensor implementation with all three eigenvalues and old-approximation comparison diagnostics. Read-only. | Supersedes old x-z principal-strain extraction; retained alongside clean derivative. |
| [02_full_tensor/implementation/project/E1_full_tensor_clean.py](02_full_tensor/implementation/project/E1_full_tensor_clean.py) | Same clean script as main reference, retained in its complete original project context. Read-only. | Full-tensor verification plus preserved clean check/smoke evidence; does not certify non-locking mechanics. |
| [02_full_tensor/implementation/project/reference_E1_original.py](02_full_tensor/implementation/project/reference_E1_original.py) | Original E1 x-z approximation baseline. Read-only. | Principal-strain extraction superseded by E1_full_tensor.py; preserve as comparison source. |
| [02_full_tensor/implementation/project/check_clean_version.py](02_full_tensor/implementation/project/check_clean_version.py) | AST preservation check with optional --smoke that executes mechanics/remodeling. Read-only development harness, not run during organization. | Compares clean/diagnostic definitions; saved smoke volume log in results/clean_smoke/. |
| [02_full_tensor/implementation/project/results/full_tensor_verification.md](02_full_tensor/implementation/project/results/full_tensor_verification.md) | Verified full-load tensor comparison, 81,531 cells. Read-only report. | Source for numerical differences in root README. |
| [02_full_tensor/implementation/project/results/principal_strain_comparison.csv](02_full_tensor/implementation/project/results/principal_strain_comparison.csv) | Raw cellwise old/full comparison. Read-only evidence. | Supports verification report; keep with phase XDMF/H5 and volume log. |

The reference and verification category folders point into the intact implementation project. Do not separate check_clean_version.py from the two scripts it reads.

## Volumetric-locking history
The following paths share prefix `03_volumetric_locking/preserved_layout/`. This complete original layout keeps sibling references and imports valid.

| Relative path under that prefix | Purpose / status (all read-only) | Later method / relevant report |
|---|---|---|
| previous_benchmark/locking_test.py | Original generated-box hydrostatic/cantilever CG1/CG2 benchmark; saved classification strong locking. | Source baseline for both mixed stages; results/locking_test_summary.md. |
| previous_benchmark/3D_ONH_E1.py | Historical original E1 constitutive reference; not the generated-box benchmark entry point. | Full-tensor reference supersedes its strain approximation; sibling ONH mesh inputs are not supplied here. |
| previous_benchmark/results/locking_test_summary.md | Original evidence and classification, hydrostatic lambda/K distinction. | Authoritative benchmark report. |
| previous_benchmark/results/locking_results.csv and hydrostatic_results.csv | Raw mechanical/hydrostatic results, with sampling and convergence columns. | Preserve plots and every raw row; never replace with new runs. |
| locking_fix/reference_E1_full_tensor.py | Immutable clean E1 snapshot, SHA-256 3b580fbb...adec530; sibling ONH mesh bundle included. | Identical to main reference; do not launch as routine inspection. |
| locking_fix/locking_plan.md | First mixed method selection and gates. | Read before changing mixed_core.py. |
| locking_fix/mixed_core.py | CG1/continuous-CG1 mixed formulation, local projection and Newton helpers. | Later interface_core.py imports helpers; continuous method fails production gate. |
| locking_fix/benchmark/locking_resistant_benchmark.py | Mixed homogeneous benchmark; reads prior CSVs. | benchmark/locking_resistant_results.csv; results/locking_fix_summary.md. |
| locking_fix/benchmark/validation_checks.py | Stabilization/tangent/energy algebra checks. | Saved validation_checks.json; not executed now. |
| locking_fix/benchmark/iterative_validation.py | Small benchmark iterative solver experiment. | Not an adopted ONH solver; saved report explains tolerances/failures. |
| locking_fix/onh_tests/ONH_CG1_baseline.py and ONH_locking_resistant.py | Mechanics-only ONH entry points. | onh_runner.py prepares reference setup; ONH_mechanics_comparison.csv and locking_fix_summary.md. |
| locking_fix/onh_tests/onh_runner.py | Reference-setup extraction, mixed mechanics and tensor/volume diagnostics. | Resolves roots before reference changes working directory. |
| locking_fix/onh_tests/ONH_mechanics_comparison.csv | Saved CG1/mixed ONH comparison. | Global improvement does not establish local accuracy. |
| locking_fix/results/locking_fix_summary.md | First mixed investigation's measured conclusions and failed production gate. | Interface investigation follows; neither is production approved. |
| locking_fix/results/CG1_steps10/ and mixed_*/ | Raw ONH fields, per-run data, volume logs and executed_reference_setup.py provenance snapshots. | Preserve all sensitivity/restart results; generated setup snapshots have historical absolute output paths and are not standalone entry points. |
| locking_fix/analyze_results.py and finalize_report.py | Historical analysis/report generators; may write outputs/reports. | Inspect as source only during documentation tasks. |
| locking_fix/interface_validation/interface_pressure_plan.md | Interface-capable pressure choice and gates. | interface_validation_summary.md records outcome. |
| locking_fix/interface_validation/interface_core.py | CG1/DG0 same-tissue jump-stabilized method; imports parent mixed_core helpers. | Pressure-jump patch passed; heterogeneous local gate failed. |
| locking_fix/interface_validation/benchmark/heterogeneous_interface_benchmark.py | Exact laminate patch and heterogeneous bending; CG1, continuous mixed, DG0 and CG2 comparisons. | heterogeneous_interface_results.csv and interface_validation_summary.md. |
| locking_fix/interface_validation/benchmark/bending_local_errors.csv | Matching-cell local J, pressure and E1/E2/E3 errors. | Decisive evidence of unresolved soft-side errors; analyze_interface_errors.py is the historical generator. |
| locking_fix/interface_validation/benchmark/homogeneous_locking_check.py and .csv | New DG0 global locking relief with saved earlier comparisons. | 81.57% nominal medium displacement-deficit recovery. |
| locking_fix/interface_validation/benchmark/stability_checks.py and .json | Pressure kernel/coupling/spectrum and tangent checks. | Limited tested stability evidence; not universal nonlinear proof. |
| locking_fix/interface_validation/benchmark/*_interface.csv, *.npz, per-case JSON | One-sided interface diagnostics, raw fields and solver histories. | Keep grouped; individual generated fields are not enumerated here. |
| locking_fix/interface_validation/interface_validation_summary.md | Authoritative failed heterogeneous-gate report and local sampling limitations. | No full ONH solution or production remodeling candidate from this method. |
| locking_fix/interface_validation/previous_stage_hashes.json | Historical relative-path preservation manifest. | Internal previous_benchmark/locking_fix layout retained to preserve its meaning. |

[Original benchmark report](03_volumetric_locking/preserved_layout/previous_benchmark/results/locking_test_summary.md), [first mixed report](03_volumetric_locking/preserved_layout/locking_fix/results/locking_fix_summary.md), [interface report](03_volumetric_locking/preserved_layout/locking_fix/interface_validation/interface_validation_summary.md).

## Coarse-CG2 direction and meshes
| Path | Purpose / status (read-only existing files) | Evidence / relationship |
|---|---|---|
| [03_volumetric_locking/cg2_tests/coarse_mesh_baseline/E1_full_tensor_clean.py](03_volumetric_locking/cg2_tests/coarse_mesh_baseline/E1_full_tensor_clean.py) | Actual coarse baseline uses CG1, not a validated CG2 candidate; differs from main reference. | Unchanged copy of archived 3D_ONH/E1_corase; sibling mesh and results retained. |
| [03_volumetric_locking/cg2_tests/README.md](03_volumetric_locking/cg2_tests/README.md) | Planned q=4/6/8 coarse-CG2 investigation; no new implementation. | Root README distinguishes saved q=6 box benchmarks from missing ONH convergence results. |
| 03_volumetric_locking/cg2_tests/coarse_mesh_baseline/results/full_tensor_clean/ | Partial saved coarse result: reference and first loading increment only in PT volume log. | No completed run/CG2 claim; preserve all XDMF/H5 files. |
| 04_meshes/original_mesh/3D_ONH_domain.xdmf + .h5; 3D_ONH_facet_region.xdmf + .h5 | Original six-tissue ONH input mesh, 81,531 tetrahedra / 15,642 nodes. | Byte copies from full-tensor project; canonical scientific geometry retained. |
| 04_meshes/coarse_mesh/3D_ONH_corase_v3.inp and domain/facet XDMF/H5 | Coarse six-tissue mesh, 18,996 tetrahedra / 3,942 nodes, 2,962 tagged facets. | XDMF metadata; source conversion_v2 bundle. Not a validated non-locking replacement. |
| 04_meshes/febio_mesh/ONH_domain.xdmf + .h5; ONH_facet_region.xdmf + .h5 | Different four-tissue FEB-derived mesh, 56,808 tetrahedra / 11,072 nodes. | Conversion report and mapping CSVs in FEB project; boundary meanings differ. |

## Mesh conversion
| Path | Purpose / status (all pre-existing files read-only) | Supersession / verification |
|---|---|---|
| [01_mesh_conversion/febio_conversion/project/README_FEB_TO_FENICS.md](01_mesh_conversion/febio_conversion/project/README_FEB_TO_FENICS.md) | FEB conversion overview and saved validation. Historical install/run commands are not current permission. | output/conversion_report.txt and .json, diagnostics/final_validation.json, element_ordering_verification.md. |
| 01_mesh_conversion/febio_conversion/project/converter/convert_feb_to_fenics.py | Converter entry point for ONH-model-remodel-test.feb; sibling modules parse, subdivide, validate and write. | 9,612 mixed source cells to 56,808 tets, original coordinates preserved; distinct from INP conversion. |
| 01_mesh_conversion/febio_conversion/project/ONH-model-remodel-test.feb | Original FEB geometry/topology/surface source. | Do not import its physics into ONH model by assumption. |
| 01_mesh_conversion/febio_conversion/project/output/ | Solver pairs, diagnostic VTUs, conversion reports, FEB/DOLFIN node and parent-element CSV mappings. | Passed saved conversion/load tests; keep mapping provenance. |
| 01_mesh_conversion/febio_conversion/project/tests/test_converter.py and test_dolfin_mesh_load.py | Historical converter/load checks; may write reports when requested. | Saved 12 unit-test and legacy DOLFIN load-test evidence; not rerun now. |
| [01_mesh_conversion/febio_conversion/project/README_FEB_SIMULATION.md](01_mesh_conversion/febio_conversion/project/README_FEB_SIMULATION.md) | Later FEB integration decisions and historical run guide. | Distinct four-tissue BC interpretation; old x-z strain trigger retained. |
| 01_mesh_conversion/febio_conversion/project/remodeling_test_3d_onh_febmesh.py | Historical full FEB integration, not current full-tensor E1 reference. | diagnostics/simulation_script.diff and metadata; simulation_output/full_run_01 is STOPPED_BY_USER, short_run_01 COMPLETED shortened schedule. |
| 01_mesh_conversion/febio_conversion/project/remodeling_test_3d_onh_febmesh_test.py | Earlier load-only integration copy, stops before scientific solve. | diagnostics/integration_copy.diff and integration_copy_test.txt. |
| 01_mesh_conversion/febio_conversion/project/remodeling_test_3d_onh_v2.py and inp_xdmf_ONH_corase.py | Preserved old model/INP converter references. | Not replaced or merged with current main model. |
| 01_mesh_conversion/inp_conversion/original_mesh/inp_xdmf_ONH_corase.py | Actual file_name=3D_ONH; converts original INP despite corase in script name. | Complete original conversion directory, including INP and mesh pairs. |
| 01_mesh_conversion/coarse_mesh_conversion/conversion_v2/inp_xdmf_ONH_corase.py | file_name=3D_ONH_corase_v3; coarse conversion source. | Sibling matching INP and written XDMF/H5 retained. |
| 01_mesh_conversion/coarse_mesh_conversion/inp_project/conversion_inp_xdmf/inp_xdmf_ONH_corase.py | Other preserved converter copy; expects coarse_v3 INP but supplied sibling INP is 3D_ONH.inp. | Pre-existing mismatch; preserve source and report, do not silently fix. |

## Additional history retained in the archive
| Path below _archive_original_projects/ | Role / status (read-only) | Distinction |
|---|---|---|
| 3D_ONH/E1/remodeling_test_3d_onh_v2.py | Older E1 branch. | Separate from later structural-weakening E1_v2 and full-tensor references. |
| 3D_ONH/E1_v2/3D_ONH_E1.py | Older x-z E1 structural model with PT volume log/results. | Full-tensor calculation supersedes approximation, while all history remains. |
| 3D_ONH/E1_v3_full tensor/ | Duplicate full-tensor project context and results. | Clean script hash equals main reference; results are retained independently, not deduplicated. |
| 3D_ONH/E3/remodeling_test_3d_onh_v2.py | Compressive E3-trigger experiment with different constitutive expressions. | No promotion to current E1 biology. |
| 3D_ONH/sig1/remodeling_test_3d_onh_v2.py | Principal tensile-stress-trigger experiment with different constitutive expressions. | No promotion to current E1 biology. |
| 3D_ONH/test1/ | Historical trial model and associated files. | Retained as separate investigation, not assumed equivalent by filename. |
| 3D_ONH/FEbio_test/ and FEbio_test_backup/ | Historical FEB mechanics variants and outputs. | Preserve both: run 20260906_225823_d31b87 records FAILED; other inspected FEbio_test summaries record COMPLETED. Consult per-run schedules/hashes. |
| 3D_ONH/*.cae, *.jnl, *.rec, *.bdf, *.inp, *.step, *.sat, *.log | Original CAD/mesh preparation and coarse v2/v3 history. | Preserved natively; not regenerated or parsed as new mesh connectivity. |
| FEBio_to_FEniCS_ONH/, full tensor upgrade/, Volume locking test/ | Complete original project snapshots. | Organized copies derive from these; original manifest verifies preservation. |

## Organization records
- [PATH_DEPENDENCIES.md](PATH_DEPENDENCIES.md): relative paths, missing input exceptions, HDF references and historical absolute paths.
- [organization_audit/original_manifest.csv](organization_audit/original_manifest.csv): pre-organization paths, sizes and SHA-256 for every original file, including existing dependencies/caches.
- [organization_audit/copy_map.csv](organization_audit/copy_map.csv): complete bundle-copy mapping.
- [organization_audit/organized_copy_manifest.csv](organization_audit/organized_copy_manifest.csv): per-file organized-to-original provenance and hash verification.
- [organization_audit/xdmf_dependencies.csv](organization_audit/xdmf_dependencies.csv): XDMF HDF file references and existence checks; no HDF dataset/FE execution test.
- [organization_audit/verification.json](organization_audit/verification.json): final original/archive and copy preservation results.
- [ORGANIZATION_REPORT.md](ORGANIZATION_REPORT.md): final tree, counts, scope and checks.


## Maintenance requirements — adopted 2026-09-10
This is a living index. Read it and README_CODEX.md before every task. Update it only for meaningful important-file additions, status changes, supersession or archival; follow the [Documentation maintenance rule](README_CODEX.md#documentation-maintenance-rule--adopted-2026-09-10).

For each affected important script, report, benchmark, reference model, production candidate, mesh, converter or verification file, record:
- Relative path and purpose.
- Current status, explicitly identifying reference, experimental, failed validation, archived or production candidate as applicable. Candidate status alone does not imply validation.
- Which file it supersedes, if any.
- Associated verification/report file, if any; state when none exists.

Retain historically important entries after supersession or archival and preserve traceability. Do not index temporary artifacts or make unnecessary updates for minor fixes, unchanged reruns or exploration without a useful conclusion. At task completion, report whether this index and README_CODEX.md were updated and summarize the changes; if neither needed updating, use the exact no-update statement in the README.

| Documentation path | Purpose | Status / relationship |
|---|---|---|
| [AGENTS.md](AGENTS.md) | Persistent workspace instructions requiring the read/assess/update/report workflow. | Active project-wide instructions; links the living-document rules. |
| [README_CODEX.md](README_CODEX.md) | Living research context, conclusions, open issues and documentation maintenance criteria. | Active project source of truth alongside relevant scripts/reports; preserves prior evidence and decisions. |
| [PROJECT_FILE_INDEX.md](PROJECT_FILE_INDEX.md) | Living navigation and important-file status/lineage. | Active index; retains historically important superseded approaches. |

## Main-model candidate added 2026-09-10
| Path | Purpose / status | Supersession / evidence |
|---|---|---|
| 05_main_model/current_candidate/coarse_cg2/E1_full_tensor_clean.py | Experimental displacement-only CG2/q4 main-model candidate, optional one-increment smoke mode, unique results directories. | Derives from user-supplied E1_corase CG1 script; does not supersede immutable references or establish production approval. TEST_REPORT.md. |
| 05_main_model/current_candidate/coarse_cg2/3D_ONH_corase_v3_domain.xdmf + .h5; 3D_ONH_corase_v3_facet_region.xdmf + .h5 | Matching copied 18,996-cell mesh bundle; original sources preserved. | source_manifest.json records source paths and SHA-256 hashes; mesh load and hash checks passed. |
| 05_main_model/current_candidate/coarse_cg2/TEST_REPORT.md; smoke_q4.log; results/cg2_q4_20260910_140713_731717/ | Initial q4 smoke PASS at 1% load, two Newton iterations; raw log and output evidence. | No full-load/remodeling or quadrature/locking validation. Retain raw results. |
| 05_main_model/current_candidate/coarse_cg2/check_preservation.py; source_manifest.json; README.md | Preservation check, provenance and candidate usage/sampling limits. | Checks passed; no prior files superseded. |
## Degree-3 performance investigation added 2026-09-10
| Path | Purpose / status | Supersession / evidence |
|---|---|---|
| 05_main_model/current_candidate/coarse_cg2/E1_full_tensor_cg2_q3.py | Experimental displacement-only CG2/q3 candidate; three early loading increments converge but warm solves remain about 91 seconds. | Latest performance trial; differs from preserved E1_full_tensor_clean.py only in default quadrature degree. Does not supersede a validated reference. Q3_PERFORMANCE_REPORT.md. |
| 05_main_model/current_candidate/coarse_cg2/benchmark_loading_cost.py; README_Q3.md | Bounded original-increment loading harness and usage; separates setup/compile, solve, diagnostics and output timings. | New performance workflow, no biological changes or remodeling execution. Q3_PERFORMANCE_REPORT.md. |
| 05_main_model/current_candidate/coarse_cg2/Q3_PERFORMANCE_REPORT.md; benchmark_q3_20260910.log; results/cg2_q3_20260910_142354_476629/ | Performance finding: approximately 96% of warm solve time in MUMPS; degree reduction did not meet practical runtime target. Raw fields, timings and source snapshots preserved. | Adds recurring-cost evidence beyond the q4 smoke report; no quadrature-accuracy validation or controlled q4 speed comparison. |
## Runtime clarification and solver guidance added 2026-09-10
| Path | Purpose / status | Supersession / evidence |
|---|---|---|
| 05_main_model/current_candidate/coarse_cg2/RUNTIME_BUDGET_AND_SOLVER_OPTIONS.md | Corrected five-hour budget, verified loaded reference BLAS/LAPACK backend and proposed solver-performance options. Read-only inspection; no new simulation or optimized solver validated. | Supersedes the one-second-target interpretation in Q3_PERFORMANCE_REPORT.md, whose numerical results remain preserved. Records runtime process-map evidence and official solver references. |
## PT-only quadratic proposal added 2026-09-10
| Path | Purpose / status | Supersession / evidence |
|---|---|---|
| 05_main_model/current_candidate/coarse_cg2/PT_QUADRATIC_FEASIBILITY.md | Proposed conforming PT-P2/rest-P1 construction, potential DOF reduction and surrounding-locking limitations. Hypothesis; no mechanics implementation/validation. | Does not supersede existing candidates. Topology evidence below; read with saved locking/interface reports. |
| 05_main_model/current_candidate/coarse_cg2/inspect_pt_enrichment.py; pt_enrichment_topology.json | Mesh-only edge-incidence count and hashed source provenance; 37,191 prospective independent vector DOFs before BCs. Executed successfully with existing DOLFIN. | Supports PT_QUADRATIC_FEASIBILITY.md; not a simulation or runtime result. |
## PT+LC regional quadratic implementation and comparison — 2026-09-10
All paths below share prefix 05_main_model/current_candidate/pt_lc_quadratic/.
| Path under prefix | Purpose / status | Supersession / evidence |
|---|---|---|
| model_pt_lc.py; model_pt_lc_buffer.py; model_driver.py | Experimental displacement-only regional-order entry points, loading comparison by default, original cycle through --full-cycle. PT+LC-only is user-selected; buffered variant deferred and preserved for later. Production unvalidated. | User selection supersedes earlier buffered preference, not the measured sensitivity in comparison_20260910/COMPARISON_REPORT.md. Immutable references remain unchanged. |
| reduced_core.py | Conforming P2/P1 constraints, reduced Newton/MUMPS equations, explicit DOF mapping; serial only. Only reduced systems factorized. | New numerical method; no mixed pressure. Small-mesh checks and two completed ONH loading comparisons below. |
| validate_reduction.py; reduction_validation.json; validation.log | Passed affine reproduction/patch, finite-difference tangent and independently assembled CG1 equivalence. | Implementation evidence; not proof of non-locking ONH mechanics. |
| validate_heterogeneous_patch.py; heterogeneous_patch_validation.json; heterogeneous_patch.log | Passed exact bonded heterogeneous laminate across the degree/material interface; displacement/J errors below 5e-12. | Additional implementation evidence, not local-accuracy convergence. |
| comparison_plan.md; README.md | Predeclared comparison screens, region definition, run options, source dependencies and limitations. | Regional workflow; preserves original biological-model meaning and explicit full-cycle scope. |
| compare_regions.py; write_report.py | Matching-cell analysis, topology/source audit, plot and report generators. | Uses completed regional results only, no unrestricted CG2 reference. |
| comparison_20260910/COMPARISON_REPORT.md; comparison.json; metrics.csv; regional_comparison.png | Completed full-load regional sensitivity comparison: local PT/LC strain differences support keeping the buffer; average volumes nearly equal. | Buffered result is not ground truth; supersedes topology-only proposal status with measured evidence. |
| results/pt_lc_q3_20260910_195133_579949/; pt_lc_loading_01.log | Raw PT+LC-only ten-increment full-loading result, snapshots/maps and cellwise fields; all increments passed. | 40,317 independent DOFs, mean solve 40.33 s per 10%-load increment; preserve raw outputs. |
| results/pt_lc_buffer_q3_20260910_195914_464727/; pt_lc_buffer_loading_01.log | Raw buffered ten-increment full-loading result, snapshots/maps and fields; all increments passed. | 45,825 independent DOFs, mean solve 54.25 s per 10%-load increment; preserve raw outputs. |

## PT+LC degree-2 test — 2026-09-10
Paths share prefix 05_main_model/current_candidate/pt_lc_quadratic/.
| Path under prefix | Purpose / status | Supersession / evidence |
|---|---|---|
| analyze_q2_test.py | Read-only analysis of completed q2/q3 states, checks source snapshots and matching topology/state metadata, writes a fresh report folder. | New analysis workflow; does not rerun mechanics or overwrite raw results. |
| q2_test_20260910_223658_922754/REPORT.md; comparison.json | q2 loading PASS; warm solve 38.30 s, saved q3 40.55 s; full-load PT/LC E1 differences 0.0161%/0.0221%, no PT threshold switches. Experimental evidence, not quadrature convergence. | Supplements q3 regional comparison; user-selected PT+LC active, buffered variant deferred. |
| results/pt_lc_q2_20260910_223658_922754/; pt_lc_q2_loading_01.log | Raw ten-increment full-loading q2 fields, timings and source snapshots, all increments passed with four Newton iterations. | Explicit --quadrature-degree 2, same 40,317 independent DOFs and unchanged biology. Retains q3 defaults/sources/results; no full remodeling test. |

## PT+LC q4 versus q2 comparison — 2026-09-10
Paths share prefix 05_main_model/current_candidate/pt_lc_quadratic/.
| Path under prefix | Purpose / status | Supersession / evidence |
|---|---|---|
| analyze_q2_vs_q4.py | Matching-state q2/q4 comparison, source-snapshot checks, local E1/J, threshold and displacement-coefficient metrics. | New analysis derived from analyze_q2_test.py; preserves earlier sources and raw results. |
| results/pt_lc_q4_20260910_224942_124198/; pt_lc_q4_loading_01.log | Raw PT+LC-only q4 loading PASS, ten increments and four Newton iterations each, source snapshots and fields. | Same regional space/biology as q2; warm mechanics 41.45 s. This is not the earlier unrestricted full-CG2 q4 smoke. |
| q2_vs_q4_20260910_223658_922754/REPORT.md; comparison.json; DECISION.md | Paired evidence and working recommendation to use explicit q2 for further experiments; q4 retained for later checks. | Full-load PT/LC E1 relative differences 0.0219%/0.0266%, no threshold switches at any tested level, q4 about 8.2% slower warm. Supplements q2/q3 evidence; remodeling accuracy/runtime remain open. |

## User-run full cycle located — 2026-09-11
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/current_candidate/pt_lc_quadratic/results/pt_lc_q2_20260910_233336_986066/ | Preserved user-run PT+LC-only q2 full-cycle outputs: loading, growth and unloading XDMF/H5 pairs, PT volume log, source snapshots and summary. | summary.json records PASS, 700 solves, 8654.28 s; volume log reaches unloading 100. First recorded completed regional full cycle; execution evidence only, field/remodeling accuracy not assessed. Supplements prior loading comparisons. |

## Three-layer cube partial-CG2 validation — 2026-09-11
Paths below share prefix 05_main_model/validation/three_layer_cube/.
| Path under prefix | Purpose / status | Supersession / evidence |
|---|---|---|
| PLAN.md; README.md | Layered cube definition, semi-analytical derivation, boundary conditions, predeclared accuracy screens, usage and failed-attempt history. | Isolated numerical validation; current ONH model/biology unchanged. |
| run_case.py; run_suite.py; reduced_core_snapshot.py | Paired analytical/bending full and partial models, two stiffness arrangements, n=3/6/9/12 refinement and n=6 q2 sensitivity. Regional core is byte-identical to active ONH core. | Native DOLFIN full-CG2 cube reference; no unrestricted full-CG2 ONH simulation. |
| results/*_v3/; suite_v3.log; *_v3.log | Four completed analytical pairs and ten bending pairs, including q2 tests; raw NPZ/XDMF/H5 fields and nonlinear histories. | Analytical consistency passes; partial/full bending accuracy screens fail. Earlier failed attempt folders/logs retained separately. |
| check_reactions.py; reaction_checks.json; reaction_check.log | Saved-field reaction audit with constraint-force transfer through P^T R; no equilibrium solve. | Authoritative reaction accounting, superseding earlier raw-ambient reaction sums in some case summaries; displacement results unchanged. |
| analyze.py; review_20260911_125815/REPORT.md; all_summaries.json; refinement.json; displacement_refinement.png | Completed review and plot: finest partial face displacement deficit 33.60% stiff-top / 5.78% soft-top; full reference scalar refinement changes <0.9%. | Adds failed practical accuracy evidence to earlier successful implementation checks. Does not directly quantify ONH error or validate local reference strains. |

## ONH regional convergence with frozen remodeled states — 2026-09-11
Paths share prefix 05_main_model/validation/onh_region_convergence/.
| Path under prefix | Purpose / status | Supersession / evidence |
|---|---|---|
| PLAN.md; README.md; study.py | Declared stability screens and completed 0-to-6-layer study at healthy and saved step-250/500 theta, identical frozen material fields across models. | Extends earlier one-layer/healthy comparison; no full-CG2 ONH or new remodeling trajectory. Existing candidate sources remain unchanged. |
| topology.py; topology.json | Read-only region-size audit: layer 7 is full CG2 and excluded; layer 6 leaves 143 linear cells. | Enforces standing full-CG2 restriction. |
| results/material_states.npz; import_audit.json; snapshot_*.py; summary.json | Frozen fields with bijective topology/tag audit, executed source snapshots, solver histories and complete metrics. | 21 equilibria completed; solver/BC/constraint checks retained. Arithmetic theta mean differs from volume-weighted log by definition. |
| results/layer_00/ through layer_06/; study_v2.log | Raw per-region/per-state NPZ and XDMF/H5 outputs, region maps and execution log. | Preserved evidence; initial missing-h5py startup log retained separately. No software installed. |
| report.py; review/REPORT.md; INTERPRETATION.md; comparison.json; region_trends.png; theta_log_audit.json | Completed regional sensitivity and practical plateau evidence; strict two-expansion gate narrowly missed by 8 PT threshold switches at 4-to-5 in final state. | Supersedes assumption that one layer adequately tests ONH regional sensitivity. Zero-to-six healthy PT u/E1 differences 10.12%/9.84%; full locking/mesh/trajectory accuracy still unverified. |
## Regional cost-versus-sensitivity figures — 2026-09-11
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/onh_region_convergence/plot_cost_tradeoff.py; cost_tradeoff/README.md; metrics.json; layer_curves.png/pdf; time_strain_tradeoff.png/pdf | Reproducible saved-result plots comparing 0–6 surrounding quadratic layers, correction time, amortized Newton cost and remaining strain sensitivity. Derived analysis only. | Supplements regional review: 3 layers conditional economical candidate, 5 more conservative for remodeling thresholds; preserves failed strict two-expansion gate and all raw results. Timings do not predict full-cycle runtime. |
## ONH v3/v4 mesh sensitivity — 2026-09-12
Paths share prefix 05_main_model/validation/mesh_v3_v4_20260912/.
| Path under prefix | Purpose / status | Supersession / evidence |
|---|---|---|
| mesh/; audit.py; audit.json; region_map.npz; PLAN.md; README.md | Preserved supplied 37,557-cell v4 input pairs, hashes, tissue/geometry audit and coarse three-layer physical-region transfer. | PT 18,235 and LC 1,243 cells; all fine centroids in matching coarse tissues. Both meshes still one LC element thick. |
| run_fine.py; run_fine.log; results/ | Completed healthy full-load q2 reduced equilibrium and source snapshots; solver checks PASS. | 135,543 independent DOFs, three Newton updates, 399.38 s correction. No full-CG2 or remodeling trajectory. Coarse saved layer-3 state remains immutable. |
| analyze.py; analyze_interior.log; write_report.py; review_interior/REPORT.md; comparison.json; common_samples.npz; strain_comparison.png/pdf | Authoritative common-point 27/64-interior-sample analysis and full-tensor evaluation checks. Local sensitivity screens FAIL: PT/LC DG0 E1 8.81%/15.07%, pointwise E1 6.54%/9.22%. | Extends same-mesh layer study: small mean differences do not establish local remodeling accuracy. Finer mesh is not ground truth; locking/thickness/trajectory convergence unresolved. |
| review/; analyze.log | Preserved preliminary 4/14-point comparison. | Superseded by review_interior/ because 14-point rule includes edge samples that overrepresent nonmatching faceted boundaries; raw equilibrium fields unchanged. |
## Partial CG2 on historical mixed-method benchmarks — 2026-09-12
Paths share prefix 05_main_model/validation/partial_previous_benchmarks/.
| Path under prefix | Purpose / status | Supersession / evidence |
|---|---|---|
| PLAN.md; README.md; run.py; run_v2.log; results_v2/ | Eight completed displacement-only regional benchmark cases, q6, theta=1; zero/three layers around left half, homogeneous nu=.49/.30 and heterogeneous nu=.49. Preserved fields, source snapshots and solver histories. | Reuses historical CG1/mixed/CG2 results read-only. No pressure added or biology changed. |
| report.py; review/REPORT.md; review/bending_comparison.png/pdf | Substantial compliance improvement; partial+3 heterogeneous displacement and both interface screens PASS. Whole stiff-region relative strain screen FAIL. | Supplements, does not supersede, earlier cube/ONH failures. Homogeneous medium deficit 3.28%, heterogeneous 0.133% versus same-mesh CG2; not general locking-free validation. |
| results/; run.log | First completed solve followed by CSV case-label lookup failure. | Preserved failed analysis attempt; corrected complete suite is results_v2/. |
## Consolidated locking/convergence resumption note — 2026-09-12
| Path | Purpose / status | Relationship |
|---|---|---|
| 05_main_model/LOCKING_CONVERGENCE_HANDOFF.md | Primary resume note: user-preferred PT+LC+3-layer/q2 economical candidate, benchmark passes and failures, unresolved ONH local accuracy, future evidence gates, runtime caveats and deferred constitutive/mixed-method ideas. | Consolidates prior reports without superseding raw evidence. Explicitly distinguishes candidate preference from unchanged zero-layer launcher and completed zero-layer trajectory. Links authoritative reports/results and preserves failed histories. |
## Standalone PT+LC+3-layer q2 development candidate
| Path | Purpose / status | Relationship |
|---|---|---|
| 05_main_model/current_candidate/pt_lc_3layer_q2/ | Self-contained coarse mesh candidate: model_pt_lc_3layer.py, unchanged reduced_core.py, four v3 mesh files, README, build_candidate.py, creation_provenance.json and setup_check.log. Default full 100/500/100 cycle, q2, three layers. | New user-requested development entry point; setup-only PASS (14,212 quadratic cells / 63,279 independent DOFs). No simulation or production validation. Preserves earlier zero/one-layer candidates; replaces their use as the starting point for the user's preferred three-layer development. |

## Passive porous cube investigation — 2026-09-14
Paths share prefix 05_main_model/validation/passive_porosity_cube/.
| Path under prefix | Purpose / status | Evidence / relationship |
|---|---|---|
| PLAN.md; run.py; report.py; finalize_completed.py | Isolated one-material porous-cube phenomenological law, analytical/3D checks and sensitivity; no Fa or ONH changes. | New constitutive hypothesis experiment; 1% seed porosity and stiffness exponents uncalibrated. Does not supersede preferred ONH candidate. |
| results/20260914_165652_557482/; run_04.log | 124 completed FE equilibria, source snapshots and XDMF/H5/NPZ fields. Original summary ends FAILED from analytical-postprocessing missing argument. | Preserve raw status/history; FE gates passed. Finalization below did not rerun FE. |
| results/20260914_165836_828952/summary.json; review_20260914_165836_828952/REPORT.md; metrics.json; comparison.png/pdf | Finalized PASS mechanism checks and scientific interpretation; 2.08%/14.52% hydrostatic volume loss at theta=.8 for two proposed laws. | Analytical agreement supports implementation only; large contraction requires extreme effective softening. Thickness/volume and unloading behavior distinguished. |
| run_01.log; run_02.log; run_03.log; earlier results directories | Preserved syntax, initial-guess/line-search and analytical-agreement failures. | Final homogeneous checks use boundary-compatible lateral lifting; no arbitrary-start or nonuniform stability claim. |

## ParaView view of porous cube — 2026-09-14
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/passive_porosity_cube/visualize_paraview.py; paraview_20260914_201702/ | ParaView 6.0 side-by-side hydrostatic r0/r1 state at time 9, actual-scale WarpByVector, original outlines, PNG, README and field checks. Derived visualization only. | cube_comparison.pvsm reads preserved raw fields; bounds/porosity checked and rendered image inspected. Supersedes preliminary paraview_20260914_201604 view/checks, whose checks were recorded before final time update. No new scientific result or solve. |

## Single-element porous test — 2026-09-14
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/passive_porosity_cube/single_element.py; export_single_element.py; single_element_results/20260914_212333_928425/ | One-tetrahedron Law B test, 13 states, source snapshots, summary, XDMF/H5 and direct deformed PVD/VTU. | PASS analytical and tangent checks; J=.854818786 at theta=.8. Reuses current cube definitions with mesh and integral normalization changes, not a new constitutive law. No ONH change. |

## Two-layer porous beam — 2026-09-14
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/two_layer_porous_beam/model_two_layer_porous_beam.py; README.md | Experimental 4x1x1 mm clamped PT/LC beam. PT initial porosity 1%, existing Law B r=1 with prescribed theta 1 to .8; LC never remodels. Default loading/remodeling/unloading and direct ParaView outputs. Built for user execution; trajectory not run. | Separate benchmark, supersedes nothing; ONH and cube references unchanged. Full CG2 small mesh/q2, follower pressure; not ONH strain-trigger evolution or revised dense-solid law. |
| 05_main_model/validation/two_layer_porous_beam/SETUP_REPORT.md; setup_checks/20260914_231704_129497/ | Passed mesh, zero-load residual, tensor mapping, tangent and reference VTU checks. Source snapshot and summary retained. | Setup-only evidence, no equilibrium solves or convergence/physical validation. Earlier failed API check and preliminary passed record preserved. |
| 05_main_model/validation/two_layer_porous_beam/USER_RUN_FAILURE.md; results/20260914_234421_388751/ | Subsequent user-run trajectory FAILED at loading step 28; 27 steps converged, no remodeling/unloading. Raw results preserved. | Updates earlier unrun status. Rejected trial cause not logged; decreasing conservative pore-gap bound suggests, but does not prove, admissibility/near-closure involvement. No diagnostic rerun or model edit. |
| 05_main_model/validation/two_layer_porous_beam/RUN_REPORT.md; results/20260914_234932_299157/; results/20260914_235102_335135/ | Requested q2 and explicit q4 execution attempts FAILED at loading steps 28 and 31 respectively. Partial ParaView outputs only. | Added failure diagnostics identify actual sampled pore-domain obstruction near clamp in q4; no remodeling/unloading results. Updates prior uncertainty without superseding preserved attempts. Constitutive/load settings unchanged; default remains q2. |
| 05_main_model/validation/two_layer_porous_beam/results/20260914_235930_156048/; review_20260915_001619_628900/; analyze_completed.py | User-selected .00125 MPa (25% original) q2 cycle COMPLETED: 700 equilibria, 74 deformed PVD frames, 993.8 s. Read-only analysis script and report/metrics verify phase/output consistency. Active script now defaults to this pressure. | Updates runnable status at reduced load; earlier high-load failures preserved. PT remodeling volume -2.724%, center vertical thickness -4.276%, interface center .03626 mm downward; unloaded geometry recovers. Existing incompressible-solid barrier unchanged. Trajectory checks passed, not physical or convergence validation. |
| 05_main_model/validation/two_layer_porous_beam/build_volumetric_variant.py; volumetric_softening/model_volumetric_softening.py; README.md; creation_provenance.json | Separate experimental volumetric-softening variant; selected progressive r_eff(theta), r=1.8/ramp power=2, evolving theta. Baseline model preserved. | Initial builder/provenance record constant-r version; later progressive selection documented in variant README/COMPARISON_REPORT.md. No isochoric, load, LC or ONH change. |
| 05_main_model/validation/two_layer_porous_beam/volumetric_softening/screens/20260915_013919_067898/; screen_endpoints.py; results/20260915_014052_965197/; results/20260915_014256_489318/ | Preserved constant-r endpoint screen (FAILED at r=1.9 after converged r<=1.8) and FAILED full paths at r=1.8 and 1.3. | Successful endpoints did not guarantee feasible theta trajectories. Superseded as runnable variants by progressive onset, retained as scientific evidence. |
| 05_main_model/validation/two_layer_porous_beam/volumetric_softening/results/20260915_014455_833701/; review_20260915_015013_367748/; COMPARISON_REPORT.md; compare_completed.py; comparison_metrics.json | Progressive variant COMPLETED 20/100/20 cycle, 302.1 s, 74 frames. PT volume loss 3.7435%, center thickness loss 6.7383%. Trajectory/output and endpoint comparison checks passed. | Greater reduction than baseline but not roughly 10%. Preserves singular pore law; no physical/mesh/quadrature validation. Reproduction command and changed continuation counts explicit in comparison report. |
| 05_main_model/validation/two_layer_porous_beam/CUBE_BEAM_COMPARISON.md | Interpretation of hydrostatic cube versus bonded/clamped beam; matched-pressure analytical cube loss 12.6801%. | Uses existing law/reports, no new FE solve. Pressure difference alone insufficient; mechanics attribution qualitative, beam convergence/locking unresolved. Supersedes no result. |
| 05_main_model/validation/passive_porosity_cube/contractility_checks/test_contractility.py; README.md; results/20260915_165341_784418/ | PASS 46 original-Law-B state checks on CG2 N2/q2 and N3/q4: unloaded perturbation recovery, held-size individual reactions, hydrostatic full cycle and separate active-stress detection control. Report, summary, source snapshots and XDMF/H5 fields. | Adds explicit contractility verification; supersedes nothing. Zero unloaded shrinkage/held-size reaction; known control detected at .002 N per face; unloading recovers geometry at theta=.8. Implementation evidence only, not biological/ONH validation. |

## Finite-solid porous transition — 2026-09-15
Paths share prefix `05_main_model/validation/two_layer_porous_beam/finite_solid_transition/`.

| Path under prefix | Purpose / status | Evidence / relationship |
|---|---|---|
| PLAN.md; README.md; material.py; model_transition.py; baseline_core.py | Separate experimental finite-compressibility and finite-pressure closure law; displacement-only beam runner, local copied helpers. | Constructed uncalibrated law; intrinsic dense bulk/shear response, no active shrinkage. Supersedes no original or ONH model. |
| test_material.py; material_checks/20260915_184254_216542/; setup_checks/finite_20260915_184142_723886/ | PASS material energy/tangent/accounting/dense-response and open/closed FE Jacobian/zero-force checks. | Implementation evidence only. Earlier material check and output-only failed result attempts preserved. |
| test_cube_transition.py; cube_checks/20260915_184846_003078/ | PASS 61-state CG2/q2 cube closure/release test; all 192 integration points close, then solid compresses. Direct deformed PVD and source snapshots. | Peak J=.763713935, Js=.964285272; unloading recovers. Adds finite-solid transition evidence, not physical validation. |
| results/finite_20260915_184425_918410/; results/singular_20260915_184439_409971/; results/finite_20260915_185052_186773/ | COMPLETED finite low-pressure, matched singular control and finite higher-pressure beam full cycles (20/100/20). | Remodeling volume losses 2.752988%, 2.742374%, 5.175767%; no fully closed beam integration points. Preserved raw fields, snapshots, histories. |
| analyze_results.py; review_20260915_215132_282290/REPORT.md; metrics.json; transition.png/pdf | Completed trajectory/output audit and comparison; finite closure works, low-pressure global contraction remains small. | No mesh/quadrature/locking or biological validation. Records changed shear bookkeeping and nonlinear-curve comparison limitations; preserves earlier evidence. |

## LC-only beam end-clamp comparison — 2026-09-16
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/two_layer_porous_beam/lc_end_clamp/model_transition.py; baseline_core.py; material.py; README.md; PLAN.md | Experimental copy of finite-solid beam with only LC end faces clamped. PLAN is preserved parent constitutive plan. | Only boundary mechanics changed; exact DOF audit and open/closed tangent checks PASS. Supersedes no prior model. |
| 05_main_model/validation/two_layer_porous_beam/lc_end_clamp/results/finite_20260915_235042_679397/ | COMPLETED 20/100/20 cycle, .00125 MPa, 322.1 s; snapshots, summary, measurements, deformed PVD. | PT remodeling volume loss 6.440118%, center thickness loss 5.015610%, unload recovery. No closed integration points; not physical/convergence validation. |
| 05_main_model/validation/two_layer_porous_beam/lc_end_clamp/compare.py; comparison.json; REPORT.md | Passed trajectory/output audit and matched boundary comparison. | Full-clamp reference volume loss 2.752988%; released PT ends materially increase contraction. Earlier results and main ONH preserved. |

## Additional bottom-face clamp — 2026-09-16
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/two_layer_porous_beam/bottom_face_clamp/model_transition.py; baseline_core.py; material.py; README.md; PLAN.md | Experimental LC-end-clamp derivative with all three bottom-face displacement components fixed. | Same material/settings; PLAN is preserved parent constitutive plan. Supersedes no model. |
| 05_main_model/validation/two_layer_porous_beam/bottom_face_clamp/results/finite_20260916_000621_917174/ | COMPLETED 20/100/20 cycle, 287.5 s; saved snapshots, fields and deformed PVD. | Boundary/Jacobian/trajectory/output checks PASS. Volume loss 4.305231%, center thickness loss 2.476369%; unload recovery. Not physical/convergence validation. |
| 05_main_model/validation/two_layer_porous_beam/bottom_face_clamp/compare.py; comparison.json; REPORT.md | Matched comparison against LC-only end-clamped beam. | Bottom fixation reduces volume loss from 6.440118%; alters lateral restraint as well as bending. No energy partition measured. Original evidence preserved. |

## Axial porous cube directional-response review — 2026-09-16
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/passive_porosity_cube/AXIAL_RESPONSE_REVIEW.md | Read-only saved-stretch and reconstructed-stress audit of original axial_r1 cube. | Confirms3.00% lateral versus3.54% axial remodeling contraction with zero lateral stress. Strong relative bulk softening implies negative unloaded effective Poisson ratio at theta=.8. Adds constitutive limitation, supersedes no raw results or prior zero-force checks; no simulation/model changes. |

## Initially 20 percent porous axial cube — 2026-09-16
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/passive_porosity_cube/initial20_axial/PLAN.md; build.py; run_initial20.py; REPORT.md; results/20260916_200200_312428/ | Completed fixed-theta axial loading test with20% initial porosity, bottom vertical support and minimal in-plane pins; direct deformed PVD. | Original r1 law retained.21 states PASS analytical/tangent/balance/reaction checks; lateral1.9384%, axial5.0575% contraction. Supports negative lateral response at fixed content, not physical validation. Supersedes no prior results. |

## Proposed kinematic FEB branch — 2026-09-20
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/kinematic_febio/PLAN.md | User-requested planning-only design for true3D FEB-mesh isotropic multiplicative atrophy, source stress-flow-theta evolution without scar/stiffening. | Pending user review; no implementation/trajectory. Corrects source growth classification to area growth, distinguishes volumetric stretch versus determinant. Preserves CMT/porosity and archives; supersedes no model. |

### Kinematic FEB plan revision — 2026-09-20
05_main_model/kinematic_febio/PLAN.md is now the accepted implementation design with E1-based evolution, not the superseded stress/flow trigger. User selected CG1, full3x3 tensor extraction, .005MPa, existing FEB BCs, isotropic cube-root Fg. Separate implementation task01a0bd91-09dd-78f3-9774-daa05179c179; awaiting implementation/verification results, no production accuracy claim. Older plan entry above records history.

### Kinematic FEB implementation — 2026-09-20
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/kinematic_febio/model_kinematic_febio.py; kinematic_core.py; README.md | Experimental true3D CG1/q4 isotropic kinematic model; explicit setup/smoke/full-cycle modes. E1 directly updates PT theta=Jg; no blood flow. | Implements accepted PLAN; supersedes no CMT/porosity candidate or reference. Full trajectory not run; locking/physical accuracy unresolved. |
| 05_main_model/kinematic_febio/mesh/; mesh_provenance.json | Local FEB mesh copies and original/copy SHA-256 evidence. | Identical to preserved 04_meshes/febio_mesh bundle;56,808 tetrahedra. |
| 05_main_model/kinematic_febio/verify_kinematic.py; verification/20260920_145652_888987/; VERIFICATION.md | Bounded analytical/FE consistency checks and report. | PASS natural shrinkage, restraint reactions, objectivity, full tensors, energy/tangent and scalar evolution. Not physical/locking validation. |
| 05_main_model/kinematic_febio/results/smoke_20260920_145833_424070/; audit_run.py | Completed10/3/10 FEB cycle, ParaView states.xdmf, source snapshots, numeric fields/logs and independent audit. | 469.89s;24 states audited PASS. Full100/500/100 cycle remains unrun. Preserves all earlier evidence. |

### Kinematic revision — 2026-09-21
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/kinematic_febio/results/full_20260920_191625_817680/ | User-run baseline full100/500/100; summary reports COMPLETED in7.79h. | Summary and loaded state inspected; entire trajectory not independently audited this task. Original sources/results immutable. Supersedes prior not-yet-run status for baseline only. |
| 05_main_model/kinematic_febio/REVISION_20260921.md; check_phase_output.py; verification/phase_output_20260921_075744_349982/report.json | Separate-phase output regression and saved-state threshold assessment PASS, no mechanics solve. | Launcher now defaults E1crit=.005, optional .02 baseline; separate phase XDMF/H5. New full trajectory untested. audit_run.py supports old/new output layouts and recorded thresholds. Preserves historical verification. |

## Confined initially porous cube — 2026-09-21
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/passive_porosity_cube/initial20_confined/run_confined.py; PLAN.md; REPORT.md; results/20260921_125051_415234/ | Completed20 loading increments with horizontal lateral constraints, unchanged20% initial porosity and original r1 law. | PASS implementation/analytical/force/tangent checks; vertical compression4.22487% versus5.05752% free-sided. Direct deformed PVD, preserved snapshots. No physical-validation claim or supersession of prior model. |

## Zero initial effective Poisson ratio cube — 2026-09-21
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/passive_porosity_cube/initial20_zero_nu/run_zero_nu.py; PLAN.md; REPORT.md; results/20260921_131755_755453/ | Completed21-state free-side axial test, G=1.5*K_initial with initial20% porosity and unchanged original volumetric law. | PASS analytical/tangent/balance/reaction checks. Axial compression11.34567%, lateral expansion1.67632%, volume loss8.34850%. Initial nu=0 only; no physical validation or prior model supersession. |

## Coupled porous-to-solid cube — 2026-09-21; rejected general material
Paths below share `05_main_model/validation/passive_porosity_cube/coupled_transition/`.

| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| material.py; run_cube.py; fe_helpers.py; PLAN.md; REPORT.md | Experimental energy-derived evolving shear/volumetric cube law, with dense-solid endpoint. **REJECTED for general use due to off-path instability.** | Open-pore zero lateral axial response demonstrated; transition calibrated to axial path and separate porosity accounting. Not an accepted beam/ONH candidate; supersedes nothing. |
| results/20260921_135909_437246/; audit_results.py; audit.json | Completed and audited 141-state fixed-porosity loading/extended-loading/unloading cube, 206.52 s. PVD/XDMF/H5 and immutable snapshots. | Same-load axial/volume loss 12.16646%, lateral approximately zero. At .025 MPa pores close, solid compresses, lateral expansion develops. Axial FE checks pass; material rejection remains. |
| verify_material.py; material_checks/20260921_140031_103231/; check_rejection.py; rejection_confirmation.json | Energy/objectivity/dense-stress checks pass; broader acoustic stability fails, confirmed independently by energy curvature. | Negative rank-one stiffness approximately -2.73575 MPa near hydrostatic J=.82742. Rules out promotion based only on attractive axial output. |
| screen_transition.py; screen_transition.csv; rejected_material_v1.py; rejected_material_v1_screen.csv; rejected_high_penalty.py; material_screen.csv | Preserved exploratory/rejected material constructions and screens. | Not alternative validated models. Startup-only failed result folders 20260921_135825_983954 and 20260921_135836_888522 are also retained. |

## Porous hyperelasticity literature — reviewed 2026-09-21
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| Literature resource/s0167-66362900064-4.pdf | Danielsson, Parks and Boyce (2004), Constitutive modeling of porous hyperelastic materials; user-supplied reference, unchanged. | Homogenized energy plus explicit-void FE comparison. Eq.31–32 are initial small-strain Mori–Tanaka moduli with incompressible matrix, not a ready finite-solid closure/remodeling law. Fig.4 tensile negative pressure differs from its descending-branch negative tangent; neither alone reproduces the rejected cube's acoustic test. Review conclusions recorded in README_CODEX.md; no implementation or simulation. |

## Mori–Tanaka cube comparisons — 2026-09-21
Paths share `05_main_model/validation/passive_porosity_cube/mori_tanaka_cube/`.

| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| run_cube.py; material.py; fe_helpers.py; verify_material.py; PLAN.md | Separate initial-modulus control and explicitly experimental current-porosity extension of Eq.31–32. | Both assume incompressible matrix for this comparison only. User's eventual finite-solid requirement retained. Not the paper's full large-strain energy; supersedes no model. |
| results/20260921_183345_512150/initial/; results/20260921_183345_512150/current/ | Two completed 101-state CG2/q4 axial loading/extended-loading/unloading cubes; raw fields and parent source snapshots. | Same-load volume losses .330845%/.310940%, lateral expansion about 1.16%. At .025 MPa volume losses 3.38142%/1.14482%. Pores remain open; no remodeling or finite-solid closure. |
| REPORT.md; audit_results.py; results/20260921_183345_512150/review/ | Passed bounded material/FE/saved-trajectory checks, comparison JSON/CSV/PNG and report. | Sampled acoustic stiffness positive. Current extension has nonmonotonic volume loss at higher axial load; physical assessment unresolved. Not global stability or biological validation. |

## Rejected coupled-transition visualization — 2026-09-21
| Path | Purpose / status | Evidence / relationship |
|---|---|---|
| 05_main_model/validation/passive_porosity_cube/coupled_transition/visualize_moduli.py; plots/20260921_202055_081852/ | Requested initial/current-porosity modulus and smooth-weight plots; PNG/SVG, CSV, source snapshot, README and checks.json. | Read-only evaluation, no FE solve. Confirms negative hydrostatic G/K, distinguishes energy coefficients from tangent moduli and two transition widths; reveals bulk tangent jumps/overshoot from full energy despite smooth quintic weights. Model remains rejected; supersedes no result. |

## Smooth saturation exploration — 2026-09-21
05_main_model/validation/passive_porosity_cube/coupled_transition/smooth_saturation/: test_curve.py, REPORT.md, checks.json, stability.csv, curves.csv, PNG/SVG. New curve/material-point trial only; smooth intrinsic endpoints achieved but finite-strain extension REJECTED by confirmed negative acoustic/energy curvature. No FE result or supersession of earlier candidates.

## Convex closure curve exploration — 2026-09-21
05_main_model/validation/passive_porosity_cube/coupled_transition/convex_closure/: test_curve.py, diagnose.py, REPORT.md, checks.json, diagnosis.json, stability.csv, curves.csv, PNG/SVG. Convex modulus proposal tested with retained energy/accounting; REJECTED near closure by independently confirmed negative energy curvature. Preserves all earlier trials; no FE trajectory.

### Kinematic storage inspection — 2026-09-22
`05_main_model/kinematic_febio/results/full_20260921_231157_207992/`: existing user-run results inspected for storage only,19.08GiB total (11.387GiB H5,7.684GiB NPZ). Every-step visualization and duplicate numeric audit arrays dominate; no deletion, settings change or new trajectory accuracy claim. See README_CODEX.md storage note. Earlier results preserved.

Storage comparison: current_candidate/pt_lc_3layer_q2/results/pt_lc_3layer_q2_20260912_224103_807304 inspected read-only:95.82MiB,72 phase frames, compact scalar/vector outputs, one NPZ; contrasts with703 frames and extensive tensor/duplicate audit output in FEB. This inspection makes no new mechanics-validation claim; see README storage comparison.

## Exploratory porous ONH preparation — 2026-09-22
05_main_model/validation/porous_onh_pt_lc/PLAN.md: user-authorized zero-layer PT+LC CG2 porous transfer with evolving theta despite failed material stability. Variant clarification pending; preparation only, no implemented solver or simulation. Existing candidates/results preserved.

### Convex porous ONH implementation — 2026-09-22
validation/porous_onh_pt_lc/ under05_main_model: model_porous_pt_lc.py, material.py, copied reduced_core.py, verify_material.py/material_verification.json, PLAN.md, README.md, VERIFICATION.md. Selected convex law with evolving PT theta, zero-layer/q2. Setup results/20260922_122731_171600 PASS; trajectory launched results/20260922_123017_015312, status in summary.json. Known negative material stiffness retained by user request; no validation or main-candidate supersession.

Porous ONH results/20260922_123017_015312/loading.xdmf (+H5) now includes reference and first converged loading frame; per-frame NPZ checkpoints and PT volume log present. Run continuing at initial observation; not a completed remodeling trajectory.

Porous ONH interruption observed2026-09-22: results/20260922_123017_015312 is incomplete, no active process; stale RUNNING summary last records remodeling188 complete/189 attempting. Latest saved checkpoint remodeling_0180.npz. No unloading; termination cause unknown. Raw files preserved; no restart during inspection.

## Porous homogenization and instability reference — reviewed 2026-09-22
Literature resource/j.jmps.2007.01.007.pdf: Lopez-Pamies and Ponte Castaneda(2007), PartI Analysis. Eq7 acoustic/strong-ellipticity criterion directly matches our constitutive diagnostic, but microstructure-based physical interpretation does not validate our constructed instability. Section2 principal-branch limitations and Section4.4 incompressible porosity identity inspected; PDF preserved. Detailed findings in README_CODEX.md. No simulation/source change.

## Porous ONH LC-bottom foundation — 2026-09-22
05_main_model/validation/porous_onh_pt_lc/lc_foundation/: model_porous_pt_lc_foundation.py, unchanged material.py/reduced_core.py copies, PLAN.md, README.md and VERIFICATION.md. User-selected normal foundation0.10MPa/mm on facet5 throughout full cycle; parent unsupported model preserved. results/20260922_170903_765060/ is SETUP_PASS only, not a trajectory. User-run command in README; no full simulation launched during preparation.

## Loaded-configuration LC foundation — 2026-09-22
05_main_model/validation/porous_onh_pt_lc/lc_loaded_foundation/: model_porous_pt_lc_loaded_foundation.py, identical parent material/reduction copies, PLAN/README/VERIFICATION. User-corrected preferred next trial: k=.10MPa/mm activates only after loading, anchored to saved loaded displacement; same anchor during unloading. Setup results/20260922_171514_404017 PASS; no full trajectory. Earlier lc_foundation (support from reference throughout loading) preserved, superseded for this request.

### Loaded-anchor cycle execution — 2026-09-22
Under05_main_model/validation/porous_onh_pt_lc/lc_loaded_foundation/: launch_detached.py; launches/20260922_172151_322486/launch.json and simulation.log; results/20260922_172152_320235/. User-authorized full cycle launched after resource/process preflight; PID2971606 confirmed active, no completion claim. Source setup previously passed; known material limitations retained.

Loaded-anchor results/20260922_172152_320235: completion confirmed from summary/log, full100/500/100 cycle, runtime5h21m13s; loading/remodeling/unloading XDMF outputs available. Process ended. Exploratory result, not independently audited physical/stability validation; earlier launched-only status superseded.

Loaded-anchor completed run summary.json and PT_volume_log.txt: endpoint comparison confirms4.54536% PT geometric volume loss during remodeling, from.7429563 to.7091863mm3. Separate from4.77685% loss relative to original reference. See README endpoint note; no new simulation.


## All-hex central rebuild — 2026-09-23
01_mesh_conversion/all_hex_core_20260923/: build.py, audit.py, README.md, inspection.json, audit.json, candidate.npz, ONH_all_hex_core.feb, ONH_all_hex_core.vtu, ONH_all_hex_boundaries.vtu and core_comparison.svg. Experimental geometry candidate: all10,476 cells hex8; unchanged9,324 outer hexes and original nodes. Topology/preservation and whole-cell positive determinant bounds pass; skewed core corners remain (minimum scaled Jacobian0.08639). Mesh-only FEB/ParaView output, not an integrated FEniCS mesh or validated simulation. Supersedes no original mesh or model.


## Whole-eye PT domain separation — 2026-09-23
01_mesh_conversion/whole_eye_pt_split_20260923/: Eye_model112_PT_separated.feb, matching visualization VTU, split_pt.py, verification.json, README.md and PT node/element correspondence CSVs. Exact original-ONH membership transfer creates separate PT domains/material (4,212 elements), leaving5,076 retinal elements. Geometry, all original element IDs/connectivity and whole-eye physics settings preserved; material clone retains source retina values. Saved-output verification PASS; no FEBio import/solve or FEniCS integration. Preserves original whole-eye/ONH sources and supersedes no baseline.


## Whole-eye smaller PT separation — 2026-09-23
01_mesh_conversion/whole_eye_small_pt_split_20260923/: Eye_model112_small_PT_separated.feb and VTU, split_pt.py, node/element correspondence CSVs, verification.json, comparison.json, README.md. New small-PT reference yields2,592 PT/6,696 remaining retinal elements, exact coordinate/connectivity matching. Preservation and VTU checks PASS; no solver/GUI verification. Separate alternative to the earlier larger-PT selection; original and previous files preserved.


## Whole-eye smaller-PT all-hex solids — 2026-09-23
01_mesh_conversion/whole_eye_small_pt_all_hex_20260923/: Eye_model112_small_PT_all_hex.feb, mesh/boundary VTUs, core_layout.svg, build.py, audit.py, README.md, audit.json, build_report.json, construction.json, quality.npz and node/core-element provenance CSVs. Geometric audit PASS:36,864hex solids, no wedges;3,060BM shells preserved. Multi-block central rebuild preserves original outer hexes, small-PT selection and settings; LC/sclera axes explicitly transferred with documented centroid sampling on new cores. Known source648 duplicate BM shell connectivities retained. No mechanics/GUI/FEniCS validation. trial_initial_block_shape/ preserved and superseded by root candidate only; earlier mixed/all-hex reference files unchanged.


### BM/fibre read-only inspection — 2026-09-23
01_mesh_conversion/whole_eye_small_pt_all_hex_20260923/shell_fibre_inspection/: README.md, inspection.json, duplicate_shell_pairs.csv, BM_locations.vtu and BM_duplicates_only.vtu. Locates648 pre-existing coincident BM pairs and explains source local fibre axes. VTU readback PASS; source/candidate unchanged, no deletion or solve. Supersedes no mesh.

## Hex-preserving whole-eye XDMF/H5 — 2026-09-23
`01_mesh_conversion/whole_eye_hex_xdmf_20260923/`: convert_hex.py, test_converter.py, verify_dolfin.py and README.md. `mesh/` contains whole_eye_domain and whole_eye_facet_region XDMF/H5 pairs; separate BM raw-shell and unique-facet pairs; mesh_metadata.npz, correspondence CSVs, conversion_report.json and MESH_HANDOFF.md. 36,864 hex8 cells preserved, PT9/LC8 tissue tags; no tetrahedralization. Format/topology/geometry roundtrip passed; five converter tests passed. WSL startup blocked the installed-DOLFIN check, so runtime/mechanics compatibility remains unverified. BM shell physics and fibre constitutive integration are not implemented by export. This is the hex-preserving alternative to the older tetrahedral conversion, not a replacement of preserved sources/results.

## BM-free whole-eye hex bundle — 2026-09-23
`01_mesh_conversion/whole_eye_hex_no_bm_20260923/mesh/`: user-selected downstream domain/facet XDMF/H5 pairs, filtered mesh_metadata.npz, correspondence CSVs, conversion_report.json and MESH_HANDOFF.md. BM shell/facet files and metadata excluded; all 36,864 hex solids and boundary data preserved exactly. FEB and older BM-containing export unchanged. No shell mechanics; runtime compatibility remains unverified. Supersedes BM-containing bundle only for the user's current integration choice.

## Whole-eye hex kinematic adapter — 2026-09-23
| Path | Purpose/status | Evidence/relationship |
|---|---|---|
| 05_main_model/kinematic_whole_eye_hex/model_whole_eye_hex.py; hex_geometry.py; README.md | Separate experimental Q1 hex growth/remodeling adapter; PT9 only, corrected z-only symmetry2, clamp1/5, IOP3/CSFP4. Compact phase outputs. | Reuses kinematic energy/evolution with correct hex integration and source mapping; previous tetrahedral/CMT/porosity candidates preserved. Runtime blocked by WSL startup; not validated or run. |
| 05_main_model/kinematic_whole_eye_hex/check_offline.py; verification/20260923_230842_775545/report.json | Offline source-hash, distorted-hex geometry, volume, affine and PT-only update checks PASS. | No DOLFIN loading/assembly/BC/solver or trajectory verification. Four material placeholders explicit in README. Uses preserved whole_eye_hex_no_bm_20260923/mesh read-only. |

2026-09-24 runtime status: README_CODEX.md records read-only WSL diagnosis (missing expected shared runtime VHDs; registered Ubuntu disks still present). No research-file or environment repair; whole-eye hex runtime verification remains blocked.

2026-09-24 WSL repair outcome: same-version signed WSL runtime repair succeeded (standard UAC, exit0); missing shared VHD files restored; Ubuntu-24.04 and legacy DOLFIN/NumPy/MUMPS import checked successfully. Logs/installer: C:/Users/E1521285/AppData/Local/Temp/codex-wsl-repair-20260924/. Supersedes startup-blocked status only; whole-eye hex model runtime/mechanics verification remains pending. See README repair record.

### Whole-eye hex runtime failure — 2026-09-24
05_main_model/kinematic_whole_eye_hex/RUNTIME_CHECK_20260924.md; probe_hex_runtime.py; results/20260924_002013_554527/; verification/hex_probe_20260924_002212_080169/ and hex_probe_20260924_002303_753588/: restored WSL runtime check FAILED actual-mesh ordering on XDMF import and Q1 direct-permutation probe. Built-in hex control PASS. Supersedes untested-runtime status with confirmed incompatibility; no mechanical trajectory. Prior geometry evidence and sources unchanged.

## Successful legacy hex ordering — 2026-09-24
| Path | Purpose/status | Evidence/relationship |
|---|---|---|
| 05_main_model/kinematic_whole_eye_hex/analyze_ordering.py; ordering_candidate.npz; audit_ordering.py | Consistent global vertex numbering with cube-symmetry and independent edge/face preservation checks. | PASS all36,864hexes; original geometry/connectivity preserved. No solver migration. |
| 05_main_model/kinematic_whole_eye_hex/verify_ordered_mesh.py; reordered_mesh/ | Verified separate domain/facet XDMF/H5 runtime bundle and hash/audit manifests. | Legacy DOLFIN Q1/label/volume/roundtrip PASS; now used by main launcher. Original input files unchanged. |
| 05_main_model/kinematic_whole_eye_hex/ORDERING_FIX_20260924.md; results/20260924_091211_449150/; check_ordering_integration.py; verification/ordering_integration_20260924_091406_813610/ | Successful model setup, BC, full tensor, boundary normal, reference-load tangent and compact-output checks. | Resolves earlier ordering blocker; no equilibrium cycle or physical/locking validation. Earlier failure report preserved as history. |
## Whole-eye hex gated execution — 2026-09-24
`05_main_model/kinematic_whole_eye_hex/run_after_smoke.py` and `launch_full.py`: persistent controller and detached full-cycle launcher, gated on completed10/3/10 short-cycle residual/determinant checks. Active short run: `results/20260924_092237_612001/` (first loading step converged; cycle incomplete). Controller: `launches/gated_20260924_092550_247714/controller.json`; full launch metadata/log path will be recorded there if the gate passes. User authorized100/500/100 full run; queued at documentation time. No completed-trajectory validation implied. Earlier results preserved.

2026-09-24 user hold: automatic whole-eye hex full-run launch CANCELLED before launch; controller PID108318 stopped and controller.json marked CANCELLED_BY_USER_NO_AUTO_LAUNCH. Short test remains running. User plans PC restart before full simulation; do not launch full until requested again.


2026-09-24: whole-eye short run results/20260924_092237_612001 COMPLETED_EXPLORATORY_NOT_VALIDATED (10/3/10,109.65min); saved residual/positive-sampled-J checks pass. Full run remains held. README_CODEX.md records read-only timing/DOF comparison to kinematic_febio full_20260921_231157_207992; timing breakdown not profiled.


2026-09-24: README_CODEX.md records installed hex quadrature counts and read-only comparison of source FEB BFGS/loading controls with the current adapter. No new numerical candidate or simulation; profiling remains pending.


2026-09-24: new05_main_model/kinematic_whole_eye_hex/q3_bfgs/ contains separate q3/BFGS model, bfgs_solver.py, check_algebra.py/algebra_report.json and README. q3 check-only results/20260924_143133_587929 PASS; algebra and single-factorization reuse verified. No mechanical cycle or speedup verified; original q4 candidate/results preserved, full run held.


2026-09-24: q3_bfgs/results/20260924_143754_479346 is the authorized matched short-cycle run (running at record time, first step91.612s). q3_bfgs/compare_runs.py and compare_when_finished.py generate saved timing/field comparison to baseline on successful completion; comparison_status.json records status. Full run held.


2026-09-24: q3_bfgs/results/20260924_143754_479346 completed; COMPARISON.md/comparison_to_q4.json report3.251x speedup and small PT/LC field differences versus q4. comparison_status.json recovered to COMPARISON_COMPLETED after native-HDF5 reader fallback; original comparison.log failure preserved. No full run launched.


2026-09-24: q3_febio_controls/ is third candidate, with bfgs_solver.py/check_controls.py and explicit README limitations. Active short run results/20260924_152339_060867; compare_three_when_finished.py/compare_runs.py generate three-setting comparison. First step passes new criteria but fails earlier strict residual diagnostic; completed comparison pending. Both baselines preserved, full run held.


2026-09-24: q3_febio_controls/results/20260924_152339_060867 completed and THREE_SETTING_COMPARISON.md plus pairwise JSON/MD reviewed. Third candidate31m12s,7.50% faster than strict q3; new convergence criteria pass,20states fail prior residual diagnostic but measured short-cycle field differences remain small. README records detailed findings; full run held.


2026-09-24: selected whole-eye candidate is q3_bfgs (second setting). Full100/500/100 results/20260924_162255_916541 launched, not completed; persistent simulation.log/launch.json in launches/20260924_162255_668153, launcher launch_full.py. Earlier hold lifted by user for this run. Third candidate remains preserved for reference.


2026-09-25: q3_bfgs/results/20260924_162255_916541 full100/500/100 COMPLETED_EXPLORATORY_NOT_VALIDATED in16h2m14s. summary.json and persistent simulation.log confirm completion; all saved residuals pass and sampled determinants positive. README records status review and limitations.


2026-09-25: 05_main_model/kinematic_whole_eye_hex/q3_bfgs/CMT_POROUS_HEX_HANDOFF.md is the transfer guide for a separate CMT porous whole-eye candidate: hex import/order, tags, BC/loading, numerical methods, constitutive replacement points, baseline evidence and verification gates. No porous whole-eye implementation or run yet.



## Porous CMT whole-eye hex transfer � 2026-09-25
05_main_model/porous_whole_eye_hex/: model_porous_whole_eye_hex.py, porous_adapter.py, unchanged material.py, hex_geometry.py, adapted strict bfgs_solver.py, verify_adapter.py, launch_smoke.py, README.md and VERIFICATION.md. Separate same-mesh/load Q1/q3 candidate; PT9 porous content law replaces Fg, E1crit.02 retained from porous source. Known unstable law unchanged by explicit selection. Setup results/20260925_141136_811705 PASS_IMPLEMENTATION; authorized10/10/10 results/20260925_141325_675732 launched, persistent launches/20260925_141325_507835/simulation.log. Preserved setup failures and all previous candidates; no full-run or physical-validation claim.

porous_whole_eye_hex/audit_run.py and audit_when_finished.py: one-shot saved-trajectory audit attached to authorized10/10/10 run; completion_status.json records pending/pass/failure. First two loading increments passed, cycle still ongoing at handoff. No full-run launch.

2026-09-25: porous_whole_eye_hex/results/20260925_141325_675732 completed10/10/10 in55m18s; trajectory_audit.json passed saved diagnostic checks, completion_status.json AUDIT_PASSED. README and VERIFICATION record phase timing and endpoint values. Supersedes running status only; known constitutive instability remains, no full run or physical validation.

05_main_model/porous_whole_eye_hex/theta_nu_comparison/: plot.py, theta_nu.png/svg and theta_nu.csv; requested J=1 comparison of implemented and proposed effective Poisson ratios versus theta. Algebraic illustration only; proposal not adopted. See README interpretation.

05_main_model/porous_whole_eye_hex/exponential_nu_candidate/: plot.py, README.md, curves.csv and exponential_theta_nu.png/svg. New requested smooth actual-porosity nu curve, compared with previous quadratic proposal; not implemented in solver. Algebraic shape and endpoint checks only.


## Exponential-nu0.1%-porosity cube - 2026-09-25
05_main_model/validation/passive_porosity_cube/exponential_nu_001/: material.py, run_cube.py, copied fe_helpers.py, check_material.py/material_checks.json, audit.py, README.md, REPORT.md; results/20260925_174413_786915 completed10/10/10, direct deformed PVD, phase XDMF, snapshots, audit.json and history.png. New E1crit.005 evolving-theta axial experiment passes bounded implementation checks but general near-closure instability persists with energy-curvature confirmation. Previous cube/ONH candidates preserved; no main-model change.


05_main_model/validation/passive_porosity_cube/exponential_nu_001/fixed_theta_08/: separate fixedtheta.8 runner and material/helper copies, README/REPORT/audit.py; results/20260925_184157_203597 completed bounded axial cycle with PVD/XDMF, NPZ and audit.json/history.png. Earlier184022 result preserved for failed final recovery check; tighter numerical residual rerun passed. Same law and loads, no ONH changes; general instability remains.


## Exponential porous whole-eye candidate - 2026-09-25
05_main_model/porous_whole_eye_hex/exponential_nu_001/: model_porous_whole_eye_hex.py, exact cube material.py, porous_adapter.py, preserved numerical helpers, verify_adapter.py, transfer_provenance.json, README.md and VERIFICATION.md. User-selected exponential law with initialporosity.001,E1crit.005,evolving PTtheta. Setup results/20260925_185444_159605 PASS, no equilibrium cycle. Old whole-eye law/results and all cube references preserved; known material instability remains.

2026-09-25: porous_whole_eye_hex/exponential_nu_001/launch_full.py launched user-authorized100/500/100 results/20260925_212858_973530 after source-hash/process/memory preflight. PID319929 active; launches/20260925_212858_725942/launch.json and simulation.log record launch/output. Not completed; previous results preserved, known instability retained.

2026-09-26: exponential_nu_001/results/20260925_212858_973530 has completed loading100/remodeling500; unloading53/100 ongoing, process active. README records end-remodeling PTvolume/theta and passing saved residuals. No full-cycle completion or validation claim.

2026-09-26: exponential_nu_001/results/20260925_212858_973530 full100/500/100 completed in18h57m17s; summary/log and ended process confirmed. All saved residuals pass, sampledJ positive. README/VERIFICATION record endpoints and unloaded geometric recovery. Supersedes ongoing status; known scientific limitations unchanged.

05_main_model/porous_whole_eye_hex/exponential_nu_001/EQUATIONS.md: source-checked constitutive/evolution/equilibrium equation reference for completed run20260925_212858_973530. Clarifies coupling and limitations; no model supersession or validation claim.

2026-09-27: README_CODEX.md records source-based clarifications for exponential_nu_001: historical shear-floor range, distinction between theta and current porosity, absent whole-eye pathwise stability audit, and untested zero-seed interpretation. Existing scripts/results and EQUATIONS.md remain unchanged.

05_main_model/porous_whole_eye_hex/exponential_nu_001/shear_curve/: plot.py, shear_vs_porosity.png/svg and curve.csv document the unchanged current shear-target curve over0-60% porosity. Algebraic plot checks passed; no mechanics or stability test.

05_main_model/validation/passive_porosity_cube/no_smoothing_exponential/: separate unsmoothed accounting cube material/runner, source provenance, README/REPORT, derivative checks, audit.py, comparison.json and two plots. Results fixed08/20260927_230048_289387 and evolving/20260927_230107_801700 include ParaView PVD/XDMF and state snapshots; both cycles completed. Baseline whole-eye retained; closure stress jump and off-path instability documented, not a promoted model.

05_main_model/validation/passive_porosity_cube/pick_freeze_exponential/: frozen-modulus experiment material.py/run_cube.py, README/REPORT, provenance, derivative checks, check_feedback.py/feedback_check.json, audit.py/comparison.json and comparison.png. Four result folders: completed fixed08_n10/evolving_n10/evolving_n40 and FAILED preserved fixed08_n40. PVD/XDMF and frozen-versus-final-target diagnostics saved. Separate investigation, not baseline replacement or validated method.

05_main_model/validation/passive_porosity_cube/zero_seed_closure/: README/REPORT, build_variants.py and comparison.json. literal_jc_s has preserved failed setup20260928_003121_054294(no equilibrium/output trajectory); control_jc_0997s has derivative checks, audit and completed fixed08/20260928_003148_581695 and evolving/20260928_003207_806364 cubes with PVD/XDMF. Baseline unchanged; exact closure plus smooth zero-seed redesign remains unimplemented.
