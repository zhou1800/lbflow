# LBFlow Publication Study

This directory contains the research draft, measured results and derivations
for the LBFlow solver. The maintained implementation and experiment
scripts are `../../solver/` and `../../scripts/`; see the repository-root README
for build instructions. The study contract is [SPEC.md](SPEC.md).

The sole author is Yu Zhou, Independent Researcher, [private contact omitted]. Code is
licensed under MIT. The article has not been submitted or accepted. Funding,
competing interests, public repository access and the final release citation
must be resolved before submission.

## Manuscript And Evidence

- [Main manuscript](paper/main.tex), [supplement](paper/supplement.tex) and
  [verified references](paper/references.bib).
- [Recorded measurement summaries](paper/data/measurements.json) and
  [full wetting table](paper/data/wetting.csv).
- [Numerical changes](notes/numerical-changes.md) and
  [forcing/relaxation derivation](notes/forcing-and-relaxation.md).
- [Original source hashes](source-manifest.json).

The recorded study used GCC 13.3, CMake 3.28.3 and Open MPI 4.1.6. Measurements
predate repository packaging: the numerical C++ sources are preserved, while
build paths and executable names changed. Raw VTI fields, checkpoints and
executable snapshots are generated evidence and are not stored in Git.
Receipt paths in the distributed summary are relative to an experiment result
root; they are not links to bundled raw files.

## Run Experiments

Run these commands from the repository root after building. Install plotting
requirements with `python3 -m pip install -r requirements-analysis.txt`.

```sh
python3 scripts/benchmarks.py verification --solver build/bin/LBFlow --root results/verification
python3 scripts/benchmarks.py capillary --solver build/bin/LBFlow --root results/capillary
python3 scripts/relax_capillary.py --solver build/bin/LBFlow --initial-root results/capillary --root results/capillary-relaxed
python3 scripts/layered_flow.py --solver build/bin/LBFlow --root results/layered-arithmetic --interpolation arithmetic
python3 scripts/layered_flow.py --solver build/bin/LBFlow --root results/layered-harmonic --interpolation harmonic
python3 scripts/wetting.py --solver build/bin/LBFlow --root results/wetting
```

The verification suite has 52 executions: 16 acceleration, 24 channel,
8 continuous/parallel/restarted runs and 4 invalid-input runs. Six comparisons
check uninterrupted versus parallel/restarted fields. Result directories must
be fresh; scripts preserve executable/input hashes and run logs. Concurrent
runs and long equilibration require adequate CPU and disk capacity.

The initial wetting sweep is not sufficient for every equilibrium claim.
`relax_wetting.py` selects a radius, optional angle/method and rank count and
continues valid complete checkpoints; failed or partial states are rejected.
For example:

```sh
python3 scripts/relax_wetting.py --solver build/bin/LBFlow --initial-root results/wetting --root results/wetting-R24 --radius 24 --ranks 4
python3 scripts/relax_wetting.py --solver build/bin/LBFlow --initial-root results/wetting --root results/wetting-R12-theta60 --radius 12 --angle 60 --end-time 48000 --ranks 2
```

Use the durations and late-change tests in the supplement when comparing to
recorded results. R12/beta=0.7/theta=60 was extended to 48000 steps; its duration
must be supplied with `--end-time 48000`. Larger continued radii use 96000,
216000 and 384000 steps for R24, R36 and R48. These are size/interface-width
sensitivity experiments at fixed lattice surface tension, not a formal
matched-physics spatial convergence study.

## Rebuild Figures And PDFs

The four figures can be regenerated directly from the included measurements,
without access to the original machine or simulation outputs:

```sh
python3 scripts/make_figures.py --paper results/regenerated-paper
```

To rebuild PDFs using Tectonic, run from `research/drlab-publication/paper/`:

```sh
tectonic main.tex
tectonic supplement.tex
```

The figure PDFs are included for LaTeX compilation. The local LaTeX sources, bibliography, figure PDFs and data form the reproducible manuscript source.
