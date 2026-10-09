# Major Revision Reproduction

Frozen experiments use source, executable and input checksums to identify
their numerical kernels. Archived results must not acquire the checksum of
a later executable.

## Reanalyze The Archived Evidence

Extract `lbflow-cpc-reproducibility.zip` into a fresh directory. Its source
tree, current PDFs, `snapshot.json`, `raw/major-revision/` and
`historical/pre-major-revision-reproducibility.zip` remain separate.
Do not edit receipts or rename individual campaign directories.

In the recorded environment (Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1), run
the frozen driver directly. Set the library paths or activate an environment
with those dependencies; the original temporary library directories are not
needed on another machine.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 /path/to/raw/major-revision/analysis-source/collect_major_revision.py --root /path/to/raw/major-revision --verify-manifest /path/to/raw/major-revision/raw-manifest.json --output reanalyzed-major-revision.json
```

Compare the output with `research/drlab-publication/paper/data/major-revision.json`.
The recorded relocated reanalysis is byte-identical. The driver checks all
3199 manifest entries, native geometry/input/executable/seed/donor hashes,
frozen reference sources and raw snapshots. It recomputes wave/window fits,
front observables, continuous-source mass budgets, manufactured errors/orders,
full-force fields and the independently constructed reduced predictions.
It does not modify the frozen input tree or rewrite receipt paths on disk.

Use `major-revision-validation.json` for the exact data/raw/PDF hashes and
the interpreted scope. A later receipt-only commit can record completed
packaging checks; the archive's `snapshot.json` identifies its scientific
source commit explicitly.

## Rebuild The New Figures And PDFs

```sh
python3 scripts/make_major_revision_figures.py --data research/drlab-publication/paper/data/major-revision.json --paper results/regenerated-major-revision
```

The generator creates `data/` and `figures/` in the destination. Five PNGs,
tables and numerical macros reproduce byte-for-byte with the recorded
Matplotlib 3.11.2 environment. PDF figure timestamps can differ; scientific
content and rendered PNGs are checked separately.
Compile `main.tex` and `supplement.tex` from their paper directory with
Tectonic or a compatible LaTeX installation. Keep all included `.tex`,
bibliography, numerical macros and figure PDFs together.

## Execute Fresh Controls

These commands create fresh calculation directories. They use the current
source; use the supplied frozen launch sources and exact archived receipts
to replay a particular historical execution. Do not silently replace an
old executable or dependency with a newer one.

```sh
python3 scripts/coupled_reference_checks.py --root results/coupled-reference
python3 scripts/phase_mechanism_controls.py --root results/phase-mechanism
python3 scripts/stokes_phase_reference.py --ny 32 --end 8000 --force 1e-7 --transformed --slip 1 --transverse-relaxation 1 --fast-stokes --root results/fourth-order-H32
python3 scripts/stokes_phase_reference.py --ny 32 --end 8000 --force 1e-7 --transformed --slip 1 --transverse-relaxation 1 --transverse-localization wall --fast-stokes --root results/wall-localized-H32
python3 scripts/periodic_contact.py waves --solver build/bin/LBFlow --root results/waves-mode3 --heights 16 32 --modes 3 --end 1000
python3 scripts/periodic_contact.py waves --solver build/bin/LBFlow --root results/waves-width6 --heights 16 32 --width 6 --end 1000
```

Dense Stokes inversion remains the default; the optional fast inverse
solves the same spatial matrix. The transverse fourth-order coefficient
is supplied by the formal expansion at the tested D=1/6, without a speed
fit. A wall/bulk diagnostic uses the fixed two-unit window and rejects
overlap that would make complementary mobility negative.
Manufactured forcing verifies consistency, while selective intervention
diagnoses numerical resistance. Neither supplies a molecular wetting law.

The main article includes five new figures. The supplement includes 15
supporting figures; all 22 old/new figure assets remain available. The
older perturbation-kernel timing data do not describe current CSF performance.
No publication, upload or payment is performed by these local scripts.
