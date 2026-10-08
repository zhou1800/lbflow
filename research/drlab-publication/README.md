# LBFlow Publication Study

This directory contains the research draft, measured results and derivations
for the LBFlow solver. The maintained implementation and experiment
scripts are `../../solver/` and `../../scripts/`; see the repository-root README
for build instructions. The study contract is [SPEC.md](SPEC.md).

The sole author is Yu Zhou, Independent Researcher, [private contact omitted]. Code is
licensed under MIT. The article has not been submitted or accepted. Funding,
competing interests, public repository access and the final release citation
must be resolved before submission.

## CPC Digital-Rock Extension

The selected target is Computer Physics Communications, Computer Programs in
Physics. The primary application is pore-scale two-phase flow on digital rock.
See the numerical scope below and the licensed
[Bentheimer example](../../examples/bentheimer48/README.md).

Use a common fresh result root, for example `results/cpc`:

```sh
python3 scripts/porous_media.py bundle --solver build/bin/LBFlow --root results/cpc/bundle --jobs 2
python3 scripts/porous_media.py entry --solver build/bin/LBFlow --root results/cpc/entry --jobs 2
python3 scripts/porous_media.py entry --solver build/bin/LBFlow --root results/cpc/entry-3d --dimensions 3 --scales 1 --jobs 2
python3 scripts/porous_media.py planar --solver build/bin/LBFlow --root results/cpc/rock-planar
python3 scripts/porous_media.py pressure --solver build/bin/LBFlow --root results/cpc/pressure-report
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset examples/bentheimer48 --root results/cpc/rock-single --mode single --ranks 2
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset examples/bentheimer48 --root results/cpc/rock-two --mode two --ranks 2
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset examples/bentheimer48 --root results/cpc/rock-contracts --mode contracts
python3 scripts/prepare_bentheimer.py --output results/cpc/bentheimer48-r2 --size 48 --refinement 2
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset results/cpc/bentheimer48-r2 --root results/cpc/rock-refinement --mode refinement --ranks 8
```

`continue_rock.py` inherits a seed's inputs and continues its exact frozen
executable from a complete outgoing checkpoint. `--check-ranks 4 8` checks
one-update global fields and every outgoing population byte before continuation.
A stopped parent's exit is retained; the completed resumed trajectory supplies
its final evidence. Fresh continuous runs need no stitching.

Build unmodified LBPM v2026.04, revision
`6d686d354e5b8140841d3601e4c8c0e4e4b77e48`, separately. The measured CPU
build used C++14/MPI and zlib, and disabled CUDA, HDF5, SILO, NetCDF and the
optional timer. Its build system also selects linker/runtime dependencies.
LBPM is an external GPL-3.0 project and is not part of the MIT code.

```sh
python3 scripts/lbpm_reference.py --solver /path/to/lbpm_permeability_simulator --source-checkout /path/to/LBPM --dataset examples/bentheimer48 --refined-dataset results/cpc/bentheimer48-r2 --root results/cpc/lbpm-reference
python3 scripts/prepare_bentheimer.py --output results/cpc/bentheimer96 --size 96
python3 scripts/performance.py --solver build/bin/LBFlow --dataset results/cpc/bentheimer96 --root results/cpc/performance --cpus 0 1 2 3 4 5 6 7
python3 scripts/make_porous_figures.py --root results/cpc --paper results/regenerated-cpc
```

Run performance after the functional campaigns, with no concurrent simulation
jobs. It requires Linux and the tested Open MPI local-rank environment. Use
CPU IDs available on your machine and record the change; timing values are
workload and machine dependent.

The three new figures can also be rebuilt from included data alone:

```sh
python3 scripts/make_porous_figures.py --paper results/regenerated-cpc
```

Together with `make_figures.py`, this regenerates all seven manuscript figures.
Text numbers and the timing table are generated from measured JSON. The raw
VTI and populations are regenerated rather than bundled. Rock results are
numerical crosschecks and transient observations, not an experimental
property, REV, general steady relative-permeability curve or validated
curved-wall wetting prediction.

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

The four original figures can be regenerated directly from the included measurements,
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
