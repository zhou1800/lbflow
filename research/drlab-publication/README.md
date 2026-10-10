# LBFlow Publication Study

This directory contains the research draft, measured results and derivations
for the LBFlow solver. The maintained implementation and experiment
scripts are `../../solver/` and `../../scripts/`; see the repository-root README
for build instructions. The study contract is [SPEC.md](SPEC.md).


## Physical Research Extension

The study addresses predictive discrete interface kinetics and capillary resistance,
supported by continuum/discrete derivations and operator-controlled
experiments. Neither a new bulk model nor a measured trapping threshold is claimed.

The major revision adds 17 full LBM executions, 30 independent FV trajectories
and seven continuous-source coupled manufactured cases. Its 30 reduced-matrix
predictions do not read evolved fronts to construct their operators; 22 have
full nonlinear-reference comparisons. Twenty-six complete phase waves span
three modes and two widths. The largest reduced mean-speed error is approximately 0.0563%;
the full-FV volume-shape error reaches 4.244%, a separate approximation limit.
The theoretical selective intervention reduces the original slip1 gaps to
1.023%, 0.452% and 0.275%, without fitting its coefficient.
See `paper/data/major-revision.json`, `major-revision-validation.json` and
[the reproduction instructions](notes/major-revision-reproduction.md).

The complete pre-revision archive is preserved in the new reproduction
package. The five new figures supplement the 17 earlier figure assets;
five appear in the revised main article and 15 in the supplement.
The primary Stokes reference keeps its dense default; `--fast-stokes`
solves the identical spatial matrix. The fourth-order and wall/bulk options
are numerical diagnostics at the tested D=1/6, not a new physical phase law.

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

Together with `make_figures.py`, `make_transport_figure.py`,
`make_capillary_figure.py`, `make_wall_figure.py` and `make_wall_flux_figure.py`,
the later core/reservoir and front-resistance figure scripts regenerate
the earlier fifteen manuscript and supplementary figures; the phase-kinetics
script below adds the latest two figure pairs.
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

The core/reservoir campaign is `paper/data/core-reservoirs.json`: 209 recorded
executions distinguish alternating force response, published core operators
and guarded pressure-driven fronts. Projection beta is matched to the entire
axial bulk recurrence, independently of pressure. Published projection
already preserves the neutral stationary ratio, so the older tangent
multiplier is cosine-specific. Signed pressure errors retain the unfavorable
refinement trend. Both source-kernel manifests in `run-kernels/` reconstruct
their frozen executables byte-identically in the original toolchain.

The flat pressure extension checks a pure force-free collar and exact
component link-flux budgets. H16/H32 pressure pairs remain bidirectional
after four durations; estimated zero drift differs from closed pressure
by at most 0.0066% of the Young reference. Front mobility is 11.3--17.8%
below separately measured pure-channel conductance. This difference is
not assigned entirely to numerical error or to a new contact-line law.
The fixed-beta refinement also shrinks physical interface thickness.
These controls do not establish a throat trapping threshold or curved-rock
wetting accuracy.

```sh
python3 scripts/published_core.py calibrate --solver build/bin/LBFlow --root results/core-calibration
python3 scripts/published_core.py slits --solver build/bin/LBFlow --root results/core-slits --calibration results/core-calibration/calibration.json --heights 16 32 48 --angles 60
python3 scripts/published_core.py channels --solver build/bin/LBFlow --root results/core-channels
python3 scripts/published_core.py drops --solver build/bin/LBFlow --root results/core-drops
python3 scripts/pressure_reservoirs.py channels --solver build/bin/LBFlow --root results/pressure-channels --tau 1
python3 scripts/make_core_reservoir_figures.py --data research/drlab-publication/paper/data/core-reservoirs.json --paper results/regenerated-core
```

`pressure_reservoirs.py fronts --seeds <held-receipts> --ratios .98 1 1.02`
uses recorded closed donors; `--closed-centered-halfwidth .00025` defines
the narrower pressure pair. `hold-fronts --seeds <front-receipts>` extends
fourfold with a fresh mass-budget baseline. `--force-removed` explicitly
changes sigma to zero as a diagnostic. `parity_forcing.py` uses the supplied
historical ghost audit/executable; `staggered_diagnostics.py` retains raw
adjacent fields. The collector documents the corrected width diagnostic in
copies of receipts while preserving all original metadata and populations.

The latest `paper/data/front-mobility.json` adds 74 executions with the
unchanged pressure-core kernel. An independent zero-flux recurrence chooses
beta for fixed measured `W09/H`; H16/H32/H48 planar production controls
verify it. Closed donors, bidirectional fronts, doubled durations and
independent pure-channel conductance separate static pressure, channel
resistance and a localized excess resistance. At H16, doubling length
changes excess resistance by 0.0333%, while the relative mobility deficit
decreases from 11.26% to 5.92%. Audited capillary-force work with independent
pure-channel velocities accounts for incremental resistance within 0.023%.
This continuum diagnostic retains absolute offsets near zero flow.
At fixed measured relative width, deficits increase from 11.26% through
18.57% to 22.30%; grid convergence and a physical contact-line law are
not established. Matching width and leading scalar diffusion does not
verify a fixed effective slip length. The paper distinguishes the known
conservative Allen-Cahn phase mobility from Cahn-Hilliard chemical mobility.

Use a fresh study root and the recorded campaign names:

```sh
python3 scripts/front_mobility.py calibrate --solver build/bin/LBFlow --root results/mobility/calibrate
python3 scripts/front_mobility.py closed --solver build/bin/LBFlow --root results/mobility/closed --calibration results/mobility/calibrate/calibration.json
python3 scripts/front_mobility.py channels --solver build/bin/LBFlow --root results/mobility/channels
python3 scripts/front_mobility.py fronts --solver build/bin/LBFlow --root results/mobility/fronts-coarse --closed-root results/mobility/closed --heights 16 32
python3 scripts/front_mobility.py fronts --solver build/bin/LBFlow --root results/mobility/fronts-fine --closed-root results/mobility/closed --heights 48
python3 scripts/front_mobility.py hold --solver build/bin/LBFlow --root results/mobility/holds-coarse --front-root results/mobility/fronts-coarse --heights 16 32
python3 scripts/front_mobility.py hold --solver build/bin/LBFlow --root results/mobility/holds-fine --front-root results/mobility/fronts-fine --heights 48
python3 scripts/front_reciprocity.py --solver build/bin/LBFlow --audit-driver build/bin/NumericalAudit --root results/mobility/reciprocal-coarse --front-root results/mobility/holds-coarse --channels-root results/mobility/channels --heights 16 32
python3 scripts/front_reciprocity.py --solver build/bin/LBFlow --audit-driver build/bin/NumericalAudit --root results/mobility/reciprocal-fine --front-root results/mobility/holds-fine --channels-root results/mobility/channels --heights 48
python3 scripts/collect_front_mobility.py --root results/mobility --output results/mobility/front-mobility.json
python3 scripts/make_front_mobility_figure.py --data results/mobility/front-mobility.json --paper results/regenerated-mobility
```

The collector recomputes diagnostics and reciprocal work from raw fields
and audit logs, with inputs, geometry, executables and checkpoints checked.

The subsequent `paper/data/phase-kinetics.json` records 67 LBM executions,
ten independent Stokes/finite-volume runs and five shared-phase ablations.
Fourteen zero-force/zero-sigma wave runs test the exact phase-stencil spectrum
at fixed physical width. Differentiating the unit normal is essential:
the predicted rates agree within 0.00074%, including a second wavelength.
The measured dependence is consistent with numerical fourth-order shape
relaxation; the known conservative Allen-Cahn mapping and signed-distance
profile-preservation property are not new results.

An optional neutral Navier closure prescribes hydraulic slip separately.
Its eighteen pure profiles, stationary phase, rotation, MPI and restart
checks pass; three default configurations retain bitwise fields/populations.
The independent Stokes solve passes analytic and manufactured controls.
Its finite-volume phase reference retains a 13.72% finest-grid dynamic
discrepancy. A shared-phase ablation reduces the mean-speed difference to
0.290%; actual-force Stokes differs by 0.0619%. This supports a discrete
kinetic contribution in this geometry; the ablation is not independent
continuum validation and no physical contact-line law is established.
An H32 cross-control exchanges the two initial planar profiles: late speeds
change by less than two parts per million. Initialization therefore does
not explain the tested dynamic discrepancy.

A separately prescribed physical slip of 4 reduces the three-grid speed
gaps to 5.40%, 3.95% and 3.23%. Six further signed pure profiles verify
the analytic wall law within 6.09e-11 relatively. The finest-grid normalized
excess resistances still differ by 22.38%; larger pure conductance partly
dilutes the speed percentage. Both prescribed-slip controls, finite-grid
pure-law offsets and final-interval changes are retained. This sensitivity
experiment does not establish a converged capillary resistance or fit slip
to front speed.

Use fresh result directories. The old no-slip campaigns use the separately
pinned pressure-core kernel; optional slip uses `navier-slip-manifest.json`.
Physical `lambda=1` corresponds to lattice `wallSlipLength=1/h`.
The separate `lambda=4` control uses `wallSlipLength=4/h`.

```sh
python3 scripts/phase_wave_spectrum.py --heights 16 32 48 --modes 1 2 --output results/phase-spectrum.json
python3 scripts/periodic_contact.py waves --solver build/bin/LBFlow --root results/phase-waves --heights 16 32 48 --end 2000
python3 scripts/periodic_contact.py waves --solver build/bin/LBFlow --root results/phase-mode2 --heights 16 32 --modes 2 --end 2000
python3 scripts/periodic_contact.py fronts --solver build/bin/LBFlow --root results/neutral-slip --heights 16 32 48 --forces -0.0000001 0.0000001 --slip 1
python3 scripts/periodic_contact.py fronts --solver build/bin/LBFlow --root results/neutral-slip4 --heights 16 32 48 --forces -0.0000001 0.0000001 --slip 4
python3 scripts/periodic_contact.py pure --solver build/bin/LBFlow --root results/pure-slip4 --heights 16 32 48 --forces -0.0000001 0.0000001 --slip 4
python3 scripts/stokes_reference_checks.py --output results/stokes-checks.json
python3 scripts/stokes_phase_reference.py --ny 32 --force 1e-7 --end 8000 --transformed --slip 1 --root results/reference-H32
python3 scripts/stokes_phase_reference.py --ny 32 --force 1e-7 --end 8000 --transformed --slip 1 --phase-stencil projection --root results/ablation-H32
python3 scripts/force_stokes_audit.py --solver build/bin/LBFlow --audit-driver build/bin/NumericalAudit --root results/force-stokes --seeds <positive-front-receipts>
python3 scripts/make_phase_kinetics_figures.py --data research/drlab-publication/paper/data/phase-kinetics.json --paper results/regenerated-kinetics
```

`collect_phase_kinetics.py --root <archived-root> --output <summary>` expects
the exact campaign topology recorded in its constants. It recomputes fields,
centroids, rates and force diagnostics, checks frozen sources/executables and
hashes retained raw files. The summary records clock/cutoff corrections;
original failed and superseded receipts remain intact and excluded from
successful comparisons. The optional slip kernel is reproducible from its
fixed base commit plus patch with byte-identical executables in the recorded
toolchain. Nine CTest groups now pass. All seventeen figure pairs can be
regenerated from included measurement summaries.
Archived receipt paths are relocated in memory under the supplied root;
original receipt and input bytes remain unchanged. Reanalysis at a different
directory was checked against the canonical summary byte for byte in the
recorded environment. This permits moving the extracted raw archive without
recreating the original temporary directory names.

The continued wall study adds `paper/data/wall-consistency.json` and
`wall-study-validation.json`. It contains independent phase recurrences,
neutral fixed points, rotation/parity and adjacent-step controls, matched
bulk drops and retained nonneutral failures. The neutral correction is
experimental; nonneutral and digital-rock wall accuracy are not established.

The subsequent `paper/data/wall-flux.json` separates complete neutral phase
transport from physical contact-angle location. Conservative phase transfer
removes resolved neutral-plane flow without filtering. Tangent, halfway
normal and published scalar ghost controls retain nonneutral failures,
bulk/contour sensitivity and adjacent-step oscillations. Closed flat D2Q9
pressure improvement does not validate reservoirs or curved rock wetting.

```sh
python3 scripts/phase_wall_flux.py --output results/phase-wall-flux.json
python3 scripts/capillary_response.py --solver build/bin/LBFlow --root results/wall-flux --heights 16 32 48 --angles 60 --walls renormalized halfway ghost --capillary-models csf --phase-boundaries bounce flux
python3 scripts/continue_menisci.py --solver build/bin/LBFlow --initial-roots results/wall-flux --root results/wall-flux-hold --walls renormalized halfway ghost --angles 60 --factor 4
python3 scripts/make_wall_flux_figure.py --data research/drlab-publication/paper/data/wall-flux.json --paper results/regenerated-paper
```

Use `wall_flux_checks.py --solver build/bin/LBFlow --audit-driver
build/bin/NumericalAudit --root results/contracts` for independent transient,
virtual-gradient, pure-phase, MPI and restart checks. Eight CTest groups include
this suite. The two new source-kernel manifests pin the frozen wall-flux/
location and ghost-comparison executables to their exact pre-patch revision;
the source patches must be applied to that revision in a separate copy.

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
