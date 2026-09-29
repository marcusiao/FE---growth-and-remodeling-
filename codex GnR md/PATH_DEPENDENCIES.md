# Path dependencies and reproducibility
Organization changed locations and added documentation only. Scientific Python, meshes, saved results and original reports were not rewritten. Read-only means project policy; original file attributes were preserved.

## Intact execution groups
- Full-tensor scripts call os.chdir(dirname(realpath(__file__))) and read sibling 3D_ONH_domain.xdmf and 3D_ONH_facet_region.xdmf. Each XDMF references its sibling H5. The organized implementation/project bundle includes scripts, check_clean_version.py, mesh inputs and saved results.
- Main-model reference contains the identical clean E1 script and its four mesh files. Its results directory would be created on execution, which was not performed. Default execution runs the original remodeling schedule and may overwrite default logs; it is not an inspection command.
- Locking code uses Path(__file__).parents and sys.path imports. Interface code imports the parent mixed_core.py; its historical previous_stage_hashes.json resolves relative to the parent of locking_fix. Benchmark comparisons read sibling previous_benchmark results. Therefore the complete layout is copied under 03_volumetric_locking/preserved_layout. Do not move interface_validation or previous_benchmark individually.
- onh_runner extracts setup from the reference, which changes working directory; restart paths are resolved before execution. The stored results/*/executed_reference_setup.py files are provenance snapshots with old absolute output directories, not portable entry points.
- Coarse E1 script changes to its own directory and expects 3D_ONH_corase_v3 domain/facet pairs. Its complete E1_corase directory is copied into cg2_tests/coarse_mesh_baseline. Its actual displacement space remains CG1.
- INP converters change directory to their own file and read file_name.inp. Original conversion uses file_name=3D_ONH; conversion_v2 uses 3D_ONH_corase_v3. Complete bundles retain their matching inputs and outputs.
- FEB converter modules/test/diagnostic scripts resolve their project root from __file__; the whole feb_xdmf layout is copied. Existing .converter-deps is historical Windows dependency content, not an installation performed here. FEB solver defaults read output/ and write a timestamped simulation_output folder relative to its script.

## Pre-existing exceptions and historical paths
- The converter under coarse_mesh_conversion/inp_project/conversion_inp_xdmf expects 3D_ONH_corase_v3.inp while that original directory contains 3D_ONH.inp. This mismatch predates organization. The complete matching coarse input/converter/output group is separately available under conversion_v2; no file was substituted.
- previous_benchmark/3D_ONH_E1.py is a constitutive reference without its expected sibling ONH mesh files. The generated BoxMesh locking_test.py does not require those external inputs. Do not equate the two entry points.
- Historical Markdown image links, run JSON and some hash manifests contain C:/Users/E1521285/Desktop/Linux_shared_file/... or /mnt/c/Users/E1521285/Desktop/Linux_shared_file/... outside the umbrella. These record former locations. They were preserved verbatim and can be stale after copying.
- Original-project relative run commands must be interpreted from their original bundle root, not the umbrella root. Do not run report generators to repair links: they can overwrite historical reports/results.
- The archive preserves internal layouts, not the former external absolute locations. No symlinks or junctions were created to impersonate those paths.
- The user-provided WSL environment was not launched or tested during organization. File availability and SHA-256 checks are not a runtime reproducibility claim.

## XDMF/H5 audit
organization_audit/xdmf_dependencies.csv enumerates each XDMF HDF filename reference and checks whether the target file exists, for archive and organized copies. Pairing was checked statically from XML; HDF datasets were not opened and mesh connectivity was not modified or recomputed. Missing references, if present, are reported as pre-existing archive issues versus organized-copy issues in verification.json.

Mesh input bundles are deliberately duplicated beside scripts and under 04_meshes. This protects relative paths; it does not authorize swapping meshes between formulations. Original/fine six-tissue, coarse six-tissue and four-tissue FEB meshes have distinct provenance and boundary semantics.

See organization_audit/path_assumptions.txt for source line matches recorded before relocation. Its absolute paths intentionally describe the initial locations.

