# epw_plrn

LiF hole polaron. `scf.in`, `ph.in`, `nscf.in` and `epw1.in` build the Wannier
data once; the polaron stages reuse the files written before them, so they must
stay in the `jobconfig` order.

| stage | args | exercises |
|---|---|---|
| `epw2.in` | 3 | SCF polaron, Gaussian init, adaptive threshold |
| `epw3.in` | 3 | `restart_plrn`, `interp_Ank_plrn`, `interp_Bqu_plrn`, `cal_psir_plrn` |
| `epw4.in` | 12 | `init_plrn = 6` from `dtau_disp.plrn`, `full_diagon_plrn` |
| `epw5.in` | 3 | `scell_mat_plrn` non-diagonal supercell SCF |
| `epw6.in` | 3 | `cal_psir_plrn` in the non-diagonal supercell |

`args = 12` is the `run-epw.sh` mode that copies `dtau.plrn` to
`dtau_disp.plrn`, so `epw4.in` starts from the displacements `epw2.in` wrote.

Only what EPW prints to stdout can be tested: `testcode` runs `extract-epw.sh`
twice in this directory, once on the benchmark and once on the new output, so
anything read out of `Amp.plrn` or `psir_plrn.*` is identical by construction.

`check_plrn_files.py` closes part of that gap. `run-epw.sh` clears the stale
polaron files before each stage and runs the checker after it, so the side files
are validated in place even though they never reach the benchmark: `Amp.plrn`
and `Ank.plrn` normalisation and the `1/N_p` convention between them,
`Bmat.plrn` and `dtau.plrn` shape and ordering, the `dos.plrn` sum rules against
both, and finite, NaN-free `psir_plrn.*`.

## Tolerances

The benchmarks are one ifort run at four ranks, which the suite reproduces
exactly; CI compares against them at its own rank count. Polaron energies drift
a few meV between toolchains, since conv_thr_plrn converges the displacement,
not the energy. Every tolerance allows at least 10 meV, twice the worst drift
over ifort and gfortran at 1, 2 and 4 ranks and over the CI container.
`eelplrn` fails first: it is the smallest number, so the same few meV is a
bigger relative error.

The polaron stages are rank-consistent: on fixed Wannier data they are
bit-reproducible from 1 to 8 ranks. What moves with the rank count is the coarse
`scf`/`ph`/`nscf` chain rerun ahead of them, in the sixth digit. `epw5.in`
full-diagonalises a small nearly degenerate supercell, which turns that digit
into ~55 meV, and `epw6.in` inherits it; hence the wider `*scell` tags. A larger
supercell barely helps -- 108 cells halves the spread and costs 11 s instead of
1 s -- so treat epw5/epw6 as a smoke test; `nrpplrn` and `check_plrn_files.py`
stay exact.

`rplrn` folds into (-0.5, 0.5] so a centre at the cell edge cannot flip between
0.0001 and 0.9999. Folded values are small -- 0.06 for epw3 -- and a centre on a
lattice site would be zero, so the centre tags use absolute error only, as do
the integer tags `ampctr` and `nrpplrn`.
