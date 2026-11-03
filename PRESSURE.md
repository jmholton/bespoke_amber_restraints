# Pressure in these MD runs

What the `PRESS` that pmemd prints means in this workflow, and how far to trust it. All tests were run in
November 2026 on 1aho, with pmemd.cuda_SPFP from Amber22 and production settings: 287 K, Langevin (`gamma_ln=1`), `cut=9`,
SHAKE on H, `netfrc=0`, OPC water. The scripts and their outputs are in
`/home2/jamesh/projects/amber/bulk_solvent/water_models/pressure_test/`, and the summaries are the `*_summary.txt` files there.

## Summary

- **Restraint forces never enter the reported pressure.** On the same coordinates and velocities, pmemd prints the
  same pressure, to 0.1 bar, with the restraints on or off.
- **`ntp=4` reads high.** leap2amber writes `ntb=2,ntp=4` for every stage, and the barometer in
  `optimize_weights_runme.com` uses it too. Against the Monte Carlo barostat, which sets the volume without using the
  virial, it reads too high by about N_water kT / V: +1300 bar in a pure OPC box and about +500 bar in the 1aho cell.
  `ntp=1` or `ntp=2` with `barostat=1, taup=99999999999` keep the volume just as fixed and read the true value.
- **Restraints do change the true pressure.** Holding the protein to its reference positions lowers the true pressure
  of the 1aho cell by about 560 bar. Released, the protein pushes on its solvent within 0.1 ps.
- **The barometer (all restraints off, `ntp=4`):** it reads high just the same, by ~500-600 bar in the 1aho cell. Releasing
  the restraints removes their real effect on the pressure, but not the `ntp=4` overcount, and the overcount does not
  depend on the time step.
- **1aho:** the cell is at a true 1 bar with about 730 waters while restrained. On the `ntp=4` scale that is a
  restrained reading of about +550 bar. The default `bulk_solvent_runme.com` cell (709 waters, `ntp=4` reading 0) was
  at about -600 bar.

`ntp=4` stays the default in this repository. To get the true pressure, subtract the correction below from an
`ntp=4` reading.

## What `ntp=4` is

The Amber22 manual (sander chapter) describes `ntp=4` as "md towards a targeted volume". Its stated uses:

- **Purpose:** "not for production but for modifying the volume of the system", for example so that REMD replicas
  share one box shape.
- **Settings:** the target is set in `&ewald` with `target_a`, `target_b`, `target_c` (cell edges) and `target_n`
  (iterations to reach it, default 100).

How our runs used it:

- **No target set:** the mdin files set none of the `target_` variables, and pmemd.cuda echoes none.
- **Volume:** it stayed at 56,230.7130 A^3 to the last printed digit for all 100 ns of the w1 run, so the mode
  holds the box where it starts.
- **Printed pressure:** this is a by-product of a volume-targeting mode, and the manual says nothing about it.

The documented way to watch the pressure without changing the volume is a constant-pressure run with a coupling time
too long to act: `ntb=2, ntp=1` (or `2`), `barostat=1`, `taup=99999999999`. leap2amber has that line commented out
above its `ntp=4` line, and the 2019 runs used the `ntp=2` form. The tests below show that this form reads the true
pressure.

## Correcting an `ntp=4` reading

    P_true ~ P_ntp4 - 138.06 x T x N_water / V        (bar; T in K, V in A^3)

At 287 K this is 39,600 bar x N_water / V. For the 1aho cell (709 waters, V = 56,231 A^3) it is 500 bar, and for
the 2542-water OPC box (V = 75,735 A^3) it is 1330 bar, where 1302 +/- 55 was measured. The rule is empirical.
It fits the two systems tested. We have not looked at the pmemd source to see why the volume-targeting mode
reports pressure this way. It
is not known whether ions or other small molecules count like waters (1aho's 26 NH4+/acetate would add ~18 bar).

## The tests

### 1. Restraint forces are not in the virial

pmemd computes `PRESS = 2/(3V) x (EKCMT - VIRIAL)` from the molecular (centre-of-mass) virial. Two checks:

- Step 1 from 22 restart files (10 frames of the 100 ns w1 run, plus prod_1 and prod_2 of six scan runs). Each was run
  with restraints on and off, through `ntp=4` and `ntp=1`: the four readings were identical for every frame.
- Two frames held almost still (dt = 1e-5 ps) for 20 steps, through `ntp=4`, `ntp=1` and `ntp=2`, restrained and
  not: the same pressure at every step.

So whatever the restraints do to the pressure, they do through the configurations they hold, not as a term in the sum.

### 2. `ntp=4` against the Monte Carlo barostat, pure OPC water

2542 OPC waters, equilibrated by the MC barostat (`barostat=2`, `mcbarint=100`) at 1 bar. From its last frame:

- 0.5 ns more of MC NPT gives <V> = 76,086 +/- 54 A^3. The volume fluctuations give kappa_T = 45.5e-6/bar, as for
  water. That makes the true pressure at the start volume (75,735 A^3) +103 +/- 16 bar.
- 0.5 ns at that fixed volume, mean after 50 ps:

| path | reads (bar) | error |
|---|---|---|
| `ntb=2,ntp=4` | +1405 +/- 53 | +1302 |
| `ntb=2,ntp=1,barostat=1,taup=99999999999` | +116 +/- 14 | +13 |
| `ntb=2,ntp=2,barostat=1,taup=99999999999` | +105 +/- 16 | +2 |

With `taup=99999999999` the volume did not change at the printed precision, so `ntp=1` and `ntp=2` are as good as
constant volume. The `ntp=4` error appears only once atoms move. Frozen coordinates (test 1) read the same through
every path, but in 20 ps continuations from the same frames, with the same random seed, `ntp=1` reads ~560 bar below
`ntp=4` from the first printed step at 0.1 ps on.

The time step does not matter. In the same box, 50 ps each:

| dt | `ntp=4` | `ntp=1` | true |
|---|---|---|---|
| 0.5 fs | +1421 +/- 25 | +121 +/- 34 | +103 |
| 1 fs | +1360 +/- 38 | +79 +/- 25 | +103 |
| 2 fs | +1479 +/- 42 | +82 +/- 24 | +103 |

### 3. The restraints' own effect, 1aho cell

The w1 cell (709 waters) at the end of 100 ns, Monte Carlo barostat at 1 and 1000 bar for 1 ns each, then read where
<V> equals the crystal's volume:

| | compressibility | true pressure at the crystal's volume |
|---|---|---|
| restrained | 8.0e-6/bar | -817 +/- 44 bar |
| all restraints off | 12.4e-6/bar | -259 +/- 51 bar |

- **Restrained MC is valid:** pmemd scales the restraint reference points with the box. The restraint energy was
  313 kcal/mol at both 1 and 1000 bar, even though the volume differs by 0.8%.
- **The restraints lower the true pressure by about 560 bar.** `ntp=4` shows the same gap: -72 bar restrained, as the
  mean over 100 ns, against +376 to +658 bar unrestrained.
- **Release time course:** in 20 ps continuations of 10 frames, releasing the restraints raises the `ntp=4` reading
  from about -130 to about +500 bar by the first printed step (0.1 ps). It stays there.
  - Over the same 0.1 ps, the protein's 1-4 electrostatic energy falls by 90 kcal/mol.
  - The restrained atoms are held off the force field's local geometry. Once free, the protein pushes harder on its
    solvent at once.
- **Which pressure describes a map:** a bulk-solvent map comes from the restrained MD, so the restrained true
  pressure is the one that describes the solvent in it.

### 4. True restrained pressure of the 1aho pressure_target scan

For each scan run, the restrained MC barostat at 1 bar ran for 1 ns from its prod_2. The true pressure was then read
off with the restrained cell's 0.449 A^3/bar. R-free is for the phenix_user_solvent example model after 3 rounds of
add_waters with the `altconf` option (bulk_solvent bench, `FPART_LINEAGE.md`).

| run | waters | `ntp=4` reading | true P (bar) | R-free |
|---|---|---|---|---|
| pressure_target -500 | 698 | -448 | -956 +/- 39 | 0.1142 |
| pressure_target 0 (default) | 709 | -11 | -624 +/- 37 | 0.1096 |
| pressure_target +500 | 719 | +297 | -466 +/- 63 | 0.1082 |
| pressure_target +1500 | 753 | +1438 | +1115 +/- 39 | 0.1074 |
| pressure_target +3000 | 789 | +3074 | +2576 +/- 45 | 0.1080 |

- **Size of the overcount:** the true values sit 320-760 bar (mean ~540) below the `ntp=4` readings, close to
  N_water kT / V.
- **Effect on R-free:** it improves until the true pressure reaches about 0 bar and is flat beyond that.

## The barometer: releasing all restraints and reading the pressure

`optimize_weights_runme.com`'s barometer releases every restraint (`ntr=0`, NMR restraints dropped) and runs
`barometer_cycles` = 10,000 steps at dt = 0.5 fs (5 ps), `ntb=2,ntp=4`. It reads the mean of the 10 printed
pressures. Run exactly that way from 10 frames of the 100 ns w1 run (709 waters), and again with `ntp=1` in place
of `ntp=4`:

| | mean over 10 frames | spread (sd) of single readings |
|---|---|---|
| `ntp=4` (as in the script) | +478 bar | 230 bar |
| `ntp=1` | -145 bar | 260 bar |
| `ntp=4` - `ntp=1`, frame by frame | +623 +/- 100 bar | |

The same cell for 1 ns with every restraint off, from two of those frames:

| start | `ntp=4` | `ntp=1` | `ntp=2` | true (MC) |
|---|---|---|---|---|
| prod_50 | +460 +/- 28 | -77 +/- 34 | +65 +/- 34 | -169 +/- 43 |
| prod_100 | +461 +/- 33 | +2 +/- 27 | -20 +/- 55 | -259 +/- 51 |

So the barometer, with all restraints released, still reads ~500 bar high through `ntp=4`:

- **What the reading means:** a barometer reading of 0 means a true pressure of about -500 bar for the unrestrained
  cell. To aim for a true 0 bar, aim the `ntp=4` barometer at about +500 bar (138.06 x T x N_water / V).
- **The `ntp=1` readings:** they sit 90-260 bar above the MC values here, against +13 bar in pure water. The MC
  numbers for the unrestrained crystal are the less certain of the two: the protein and lattice can relax during
  the 1 ns runs.
- **Noise:** a single 5 ps barometer reading scatters by about +/-230 bar from frame to frame. That is ~8 waters'
  worth in this cell, so averages over several barometer runs are needed for a decision at the 100 bar level.

## Measuring the true pressure yourself

- Directly: run the fixed-volume segment with `ntb=2, ntp=2, barostat=1, taup=99999999999` (or `ntp=1`) in place of
  `ntp=4`. The mean `PRESS` is the true pressure.
- Without the virial: run the Monte Carlo barostat (`ntb=2, ntp=1, barostat=2, mcbarint=100`) at `pres0=1` and at
  `pres0=1000`, 1 ns each, and take <V> after 0.2 ns. The true pressure at the crystal's volume V0 is
  1 + 999 x (<V>_1 - V0) / (<V>_1 - <V>_1000). This works with restraints on, since the references scale with the box.
