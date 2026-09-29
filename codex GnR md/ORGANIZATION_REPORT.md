# Organization completion report
Completed 2026-09-10. This task only inspected, preserved, organized and documented existing research.

## Final folder tree
Generated fields, dependency internals and individual historical runs are intentionally collapsed.
```text
ONH_FEniCS_Project/
|-- 01_mesh_conversion/
|   |-- febio_conversion/project/          (complete original FEB project)
|   |-- inp_conversion/original_mesh/     (converter, INP, XDMF/H5)
|   |-- coarse_mesh_conversion/
|   |   |-- conversion_v2/                (matching coarse_v3 bundle)
|   |   `-- inp_project/                  (other converter copy, unchanged)
|   `-- archived_converters/README.md     (navigation to immutable bundles)
|-- 02_full_tensor/
|   |-- reference/README.md
|   |-- implementation/project/           (scripts, meshes, results together)
|   `-- verification/README.md
|-- 03_volumetric_locking/
|   |-- original_benchmark/README.md
|   |-- mixed_continuous_pressure/README.md
|   |-- interface_validation/README.md
|   |-- preserved_layout/
|   |   |-- previous_benchmark/
|   |   `-- locking_fix/
|   |       |-- benchmark/
|   |       |-- onh_tests/
|   |       |-- results/
|   |       `-- interface_validation/
|   `-- cg2_tests/
|       |-- README.md
|       `-- coarse_mesh_baseline/         (saved CG1 script, mesh, partial results)
|-- 04_meshes/
|   |-- original_mesh/
|   |-- coarse_mesh/
|   `-- febio_mesh/
|-- 05_main_model/
|   |-- reference/                        (clean E1 and required mesh files)
|   `-- current_candidate/README.md       (no approved production candidate)
|-- _archive_original_projects/
|   |-- 3D_ONH/
|   |-- FEBio_to_FEniCS_ONH/
|   |-- full tensor upgrade/
|   `-- Volume locking test/
|-- organization_audit/
|-- AGENTS.md
|-- README_CODEX.md
|-- PROJECT_FILE_INDEX.md
|-- PATH_DEPENDENCIES.md
`-- ORGANIZATION_REPORT.md
```

## Preservation and verification
- Original four directories were relocated intact into _archive_original_projects after recording pre-organization hashes. No original research file was deleted or overwritten.
- 6,128 original files were checked against the manifest; all 6,128 archived files match SHA-256. This includes 2,007 Python files, many belonging to pre-existing vendor dependencies; no scientific Python file was modified.
- 5,840 organized research/dependency copies match their original SHA-256. Seven complete source bundles have zero missing files.
- 9,288 HDF references extracted from archive/organized XDMF files resolve to existing files. No missing HDF file reference was found. These are reference occurrences, not unique H5 files.
- No software was installed. Existing dependency trees were copied as historical project content.
- No simulation, mesh conversion, scientific test harness or scientific Python script was run. No mesh, material, BC, theta, quadrature, FE space, solver or remodeling formulation was changed.
- Verification is file integrity and static dependency checking, not a new runtime or scientific validation.
- Archive contents have no new documentation injected into them and no read-only attributes were changed. Immutability is explicitly established in root/folder documentation and AGENTS.md.

Machine-readable evidence: [verification.json](organization_audit/verification.json), [original_manifest.csv](organization_audit/original_manifest.csv), [organized_copy_manifest.csv](organization_audit/organized_copy_manifest.csv), [bundle_completeness.json](organization_audit/bundle_completeness.json), [xdmf_dependencies.csv](organization_audit/xdmf_dependencies.csv). The source-to-destination directory mapping is in [copy_map.csv](organization_audit/copy_map.csv).

## Scientific status recorded, not changed
Main current E1 reference: [05_main_model/reference/E1_full_tensor_clean.py](05_main_model/reference/E1_full_tensor_clean.py). This is the identical clean full-tensor script used by the first mixed investigation, with sibling mesh pairs supplied. Full-tensor extraction is preferred; its CG1 mechanics remains locking-sensitive.

Unresolved problem: both stabilized low-order mixed candidates substantially reduce global locking but do not establish reliable local J/E1 mechanics at heterogeneous interfaces. Removing pressure continuity alone did not solve the soft-side local error. Neither is approved for production remodeling.

Current coarse-CG2 direction: 18,996 tetrahedra / 3,942 nodes, prospective explicitly controlled quadrature convergence (q=4,6,8 proposed), using CG2 first as a possible non-locking reference. The saved coarse script is CG1 and its log is partial; no coarse-ONH CG2 quadrature validation or hundreds-of-points warning log was found.

## Documentation created
Root:
- [README_CODEX.md](README_CODEX.md): concise research context, evidence, model meaning, decisions, rules and open issues.
- [PROJECT_FILE_INDEX.md](PROJECT_FILE_INDEX.md): research navigation, read-only status, supersession and report associations.
- [PATH_DEPENDENCIES.md](PATH_DEPENDENCIES.md): preserved execution layouts and pre-existing missing/historical paths.
- [AGENTS.md](AGENTS.md): directs future Codex tasks to read the root context and respect immutable references.
- [ORGANIZATION_REPORT.md](ORGANIZATION_REPORT.md): this completion record.

Folder navigation/status READMEs:
- 01_mesh_conversion/README.md
- 01_mesh_conversion/archived_converters/README.md
- 02_full_tensor/reference/README.md
- 02_full_tensor/verification/README.md
- 03_volumetric_locking/README.md
- 03_volumetric_locking/original_benchmark/README.md
- 03_volumetric_locking/mixed_continuous_pressure/README.md
- 03_volumetric_locking/interface_validation/README.md
- 03_volumetric_locking/cg2_tests/README.md
- 04_meshes/README.md
- 05_main_model/reference/README.md
- 05_main_model/current_candidate/README.md

Audit files are new organization evidence, not simulation results. The path-assumptions text records source-line matches at the original locations.

## Future task workflow
Read README_CODEX.md first. Use PROJECT_FILE_INDEX.md to locate the relevant original report and implementation. Work only in the user-specified subfolder, using the rest as read-only context. A new validated result should update the root context/index while preserving historical evidence. Do not start a numerical investigation merely because organization is complete.

