# LBFlow Publication Study

This directory contains the research draft, measured results and derivations
for the LBFlow solver. The maintained implementation and experiment
scripts are `../../solver/` and `../../scripts/`; see the repository-root README
for build instructions. The study contract is [SPEC.md](SPEC.md).

The sole author is Yu Zhou, Independent Researcher, [private contact omitted]. Code is
licensed under MIT. The article has not been submitted or accepted. Funding,
competing interests, public repository access and the final release citation
must be resolved before submission.

## Physical Research Extension

The study addresses finite-interface transport and wetting effects, supported by a
continuum momentum derivation and a new controlled campaign. See
[the physical research direction](notes/physics-direction.md) for literature
overlap, completed evidence and the remaining capillary-trapping experiments.
Neither a new bulk model nor a measured trapping threshold is claimed.

```sh
python3 scripts/interfacial_transport.py --solver build/bin/LBFlow --root results/interfacial-transport --jobs 2
python3 scripts/reference_sensitivity.py --data results/interfacial-transport/interfacial-transport.json --output results/interfacial-transport/reference-sensitivity.json
python3 scripts/make_transport_figure.py --data results/interfacial-transport/interfacial-transport.json --sensitivity results/interfacial-transport/reference-sensitivity.json --paper results/regenerated-transport
```

The campaign contains three independent planar calibrations and eighteen
layered-flow calculations. Those measurements retain the unchanged default
operators. The sensitivity postprocessing uses SciPy from the existing
analysis requirements and compares specified reference reconstructions;
it does not isolate lattice truncation error.
Rebuild its figure from the included data with:

```sh
python3 scripts/make_transport_figure.py --paper results/regenerated-transport
```

## Closed Capillary And Model Controls

The continued study compares direct bulk pressure, residual velocity and
late geometry, including an optional local tanh reconstruction and an
established curvature/Guo surface force. The options remain experimental.
They do not establish a new multiphase theory, an exact implementation of a
published complete wetting scheme, or a pore-throat trapping threshold.

```sh
python3 scripts/capillary_response.py --solver build/bin/LBFlow --root results/menisci-base --heights 16
python3 scripts/capillary_response.py --solver build/bin/LBFlow --root results/menisci-tanh --heights 16 --walls tanh
python3 scripts/capillary_response.py --solver build/bin/LBFlow --root results/menisci-csf --heights 16 --walls legacy renormalized tanh --capillary-models csf
python3 scripts/capillary_response.py --solver build/bin/LBFlow --root results/menisci-refined-linear --heights 32 48 --angles 60 90 --walls renormalized --capillary-models perturbation csf
python3 scripts/capillary_response.py --solver build/bin/LBFlow --root results/menisci-refined-legacy --heights 32 48 --angles 60 90 --walls legacy
python3 scripts/continue_menisci.py --solver build/bin/LBFlow --initial-roots results/menisci-base results/menisci-refined-legacy --root results/menisci-legacy-hold --walls legacy
python3 scripts/continue_menisci.py --solver build/bin/LBFlow --initial-roots results/menisci-base results/menisci-csf results/menisci-refined-linear --root results/menisci-linear-hold
python3 scripts/bulk_capillary_controls.py --solver build/bin/LBFlow --root results/bulk-controls
```

Heights 16/32/48 use diffusive surface-tension and duration scaling with
fixed beta. This is a joint grid/sharp-interface sequence, not fixed-interface
grid convergence. Continuations keep complete mass/field histories.
The wall-profile control sets sigma to zero to distinguish color segregation
from capillary flow. The supplied measurements retain unfavorable results
and early drop measurements, including the reason for longer relaxation.
Use `make_capillary_figure.py --paper results/regenerated-capillary` to rebuild
the figure, numerical macros and table from the distributed summary.
`collect_capillary_controls.py` reanalyzes and checksums the separate run roots.
The [run-kernel patch](run-kernels/manifest.json) reconstructs the frozen
pre-storage-refactor source used for the refined/continued controls. Its
binary was rebuilt byte-identically in the original environment, then the
current source/binary was restored. Apply that patch only in a separate copy
with `git apply --unidiff-zero` when reproducing that historical executable; binary hashes are not portable
between toolchains.

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

Together with `make_figures.py` and `make_transport_figure.py`, this regenerates
all eight manuscript figures.
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

The continued wall study adds `paper/data/wall-consistency.json` and
`wall-study-validation.json`. It contains independent phase recurrences,
neutral fixed points, rotation/parity and adjacent-step controls, matched
bulk drops and retained nonneutral failures. The neutral correction is
experimental; nonneutral and digital-rock wall accuracy are not established.

```sh
python3 scripts/phase_recurrence.py --output results/phase-recurrence.json
python3 scripts/wall_consistency.py --solver build/bin/LBFlow --root results/neutral-wall
python3 scripts/wall_consistency.py --solver build/bin/LBFlow --root results/neutral-refinement --heights 16 32 48 --betas 0.7 --methods tangential --sigmas 0.005 --curvature-stencils normal gradient
python3 scripts/make_wall_figure.py --data research/drlab-publication/paper/data/wall-consistency.json --paper results/regenerated-paper
```

`scripts/replay_recorded_case.py` relocates an archived `input.prm`, raw
geometry and outgoing seed into a fresh directory, with optional bitwise
field comparison. Select the corresponding source kernel using the manifests
and zero-context patches under `run-kernels/`; historical results do not
silently acquire the current executable's hash. Numerical audit output
identifies root-local dimensions and origin when inspecting a parallel state.

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
