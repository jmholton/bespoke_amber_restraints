# 1aho His / C-terminus / salt test scripts (2026-11-24 to 26)

Copies of the scripts in `/scratch/jamesh/projects/amber/amide/` that `WEIGHTS_AND_B.md` and
`../MD_BULK_SOLVENT_JOINT.md` cite. They hard-code that site's paths (`1aho_hurry_test`, `run_local.csh`, the Amber and
CCP4 installs) and are kept as the record of how each run was built, launched and scored, not as portable tools.

| script | what |
|---|---|
| `launch_restraint_classes.sh` | the five restraint-class arms `ff19_opc_rc_*` |
| `setup_his_builds.sh`, `build_his.csh` | the no-salt His tautomer builds `build_2s2r_amide_h*` (source of the s68 builds' protonation and tleap stubs) |
| `salt068_base.csh` | the shared 0.68 M ammonium acetate placement (11 NH4 + 15 ACY per cell) |
| `setup_s68_builds.sh`, `build_s68.csh` | the six charge-test builds `build_s68_{cEE,EE,ED,DD,MD,PD}` |
| `launch_s68_tests.sh` | the charge-test runs `ff19_opc_s68_{cEE,EE,ED,PD}` |
| `setup_s68_flip.sh` | `build_s68_PDp`: PD with His54's ring relabelled to the P orientation (positive control) |
| `launch_s68_controls.sh` | launches runs given as name:build (the controls `ff19_opc_s68_{PDp,DD,PD_r2,ED_r2}` by default) |
| `setup_s68_more.sh` | `build_s68_{PE,DE,ME}`: the rest of His54 {HIE, HID, HIP, mix} x His64 {HIE, HID} |
| `setup_s68_samesalt.sh` | `build_s68_{QD,MD,PD,PE}_u`: His54 HIP in 1, 2 or 4 copies at the SAME salt as DD (charge left to PME's uniform background); `build_s68_{DE,DD}_n`: NH4+ on crystal water HOH51 at His64 ND1 |
| `salt_base.csh`, `setup_salt_series.sh`, `add_counter_acy.py` | salt series: ammonium acetate at 0.17 / 0.34 / 1.36 M, builds `build_s{017,034,136}_{DD,DE}` |
| `worst_restraints.py` | per-atom restraint energy w d^2 and weight over iterations, pooled over copies; waters by nearest protein atom |
| `nh4_site.py` | NH4+ and nearest water distance to His64 ND1 / Tyr42 O per copy per iteration (`amber_N.rst7`); output `nh4_site_s68_n.txt` |
| `release_stats.py` | which restraints the release step seeds and resets, per run |
| `his_neighbours_md.py` | partners within 3.2 A of the His ring N over a run's saved frames |
| `his_check.py` | His54/His64 labels, ring H and charges per copy in a build's `xtal.prmtop` |
| `his_chi.py` | His chi1/chi2 (ring orientation) in PDB files |
| `his_chi_nc.py` | His chi1/chi2 over a run's saved `amber_N.nc` |
| `his_fofc.py` | iteration-averaged Fo-Fc at the His ring atoms and ring-H sites, and the local ring Fo-Fc rms (ccp4-python) |
| `gpu_nc2mtz_test.csh`, `cmp_render.py`, `avg_frames.py` | the gemmi vs GPU render comparison (`WEIGHTS_AND_B.md`, 2026-11-25) |
