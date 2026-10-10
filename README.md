# LBFlow

LBFlow is an MPI lattice Boltzmann simulator for pore-scale flow. The maintained C++11 implementation supports single-phase
flow and immiscible two-phase color-gradient flow on D2Q9 and D3Q19 lattices,
with SRT/MRT collision, voxel solids, geometric contact angles, VTI output and
population checkpoints.

Build the solver from the repository root.

## Build

Requirements: a C++11 compiler, CMake 3.16 or newer and MPI development libraries.
Tests also need Python 3.10 or newer and NumPy. On Ubuntu:

```sh
sudo apt-get install g++ cmake libopenmpi-dev openmpi-bin python3-numpy
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

The executable is `build/bin/LBFlow`. To build only the simulator without Python
or NumPy, configure with `-DBUILD_TESTING=OFF`. Builds do not download
third-party source or modify the source tree. The supported and tested platform
is Linux with MPI; other platforms have not yet been verified.

## Run A Small Volume

```sh
python3 examples/generate.py channel --output results/channel
(cd results/channel && mpiexec -n 2 ../../build/bin/LBFlow input.prm)
```

This creates and simulates an 8 x 10 x 4 voxel channel. See
[`examples/README.md`](examples/README.md) for the analytical reference and a
small two-phase drop. Run the simulator from the input directory so its relative
geometry and output paths resolve correctly. Use ParaView to view `field.*.vti`.

For installation, use `cmake --install build --prefix /your/install/directory`.

## Numerical Verification

Verification covers MRT phase viscosity and force integration,
physical velocity reporting, checkpoint time/stage semantics, field-based
convergence and numerical failure exits. The color-gradient
perturbation and recoloring remain the bulk mechanism. Optional harmonic
viscosity and wall-gradient reconstruction have separate measured tradeoffs.
An experimental curvature surface force supports closed equal-density color
flow and guarded flat pressure-reservoir experiments. Independently matched
published recoloring/MRT controls separate core and wall effects. The force
reduces residual velocities in the recorded slit controls,
with remaining pressure and wall-profile errors documented explicitly.
Fixed-width, length and rate controls trace a localized excess front
resistance to independently weighted capillary force. That response still
varies with grid resolution and is not a validated physical contact-line law.
See the [numerical change log](research/drlab-publication/notes/numerical-changes.md).

The [numerical mechanism study](research/drlab-publication/notes/major-revision-theory.md)
adds an exact 26-wave phase-spectrum comparison, a theoretical fourth-order
reference intervention and a predictive low-force capillary-resistance model.
Coupled continuous-source verification, full force/velocity fields and
wall/bulk/surface-tension controls bound this numerical mechanism. The
production C++ kernel is unchanged; no physical contact-line law is claimed.

CTest exercises coefficients from the actual C++ implementation, 16 analytic
force cases, analytic phase fluxes in isolated pores, capillary-entry direction,
single/two-rank and restarted field agreement, and invalid parameter/checkpoint
failures. A GitHub Actions workflow is provided for GCC and Clang; its execution
depends on the repository's Actions settings.
The fourth CTest group checks CSF equilibrium/forcing in both lattices and
collision operators, normal-field MPI exchange, restart and unsupported-input
failures.

The extended study has 24 analytic channel cases, 12 layered-flow cases and 30
wall-drop configurations. Its completed quantitative scope is equal reference
densities, standard pressure weights, low speed and relaxation times 0.6--1.4.
Curved-wall wetting accuracy, moving contact lines and high density ratios need
further validation. MPI is implemented; CUDA and OpenMP acceleration are future work.

## Pore-Scale Digital Rock

The [Bentheimer example](examples/bentheimer48/README.md) includes small,
checksum-pinned crops of a DOI-cited public sandstone image, with its separate
CC BY 4.0 data license and full preparation manifest. Workflows measure
single-phase permeability, phase-fraction-weighted advective Darcy flux,
component mass, late field changes and MPI/restart agreement. The study examines numerical interface kinetics and capillary resistance,
with rock flow as a supporting demonstration.

```sh
python3 scripts/porous_media.py bundle --solver build/bin/LBFlow --root results/pore-bundles
python3 scripts/porous_media.py entry --solver build/bin/LBFlow --root results/capillary-entry
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset examples/bentheimer48 --root results/rock-single --mode single --ranks 2
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset examples/bentheimer48 --root results/rock-two --mode two --ranks 2
```

The small periodically repeated image is not a representative elementary
volume. Many pores are comparable to the diffuse interface width. Transient
phase flux is reported directly and is not treated as steady relative
permeability. See the [study scope](research/drlab-publication/README.md).

## Reproduce The Study

```sh
python3 -m pip install -r requirements-analysis.txt
python3 scripts/benchmarks.py verification --solver build/bin/LBFlow --root results/verification
python3 scripts/benchmarks.py capillary --solver build/bin/LBFlow --root results/capillary
python3 scripts/wetting.py --solver build/bin/LBFlow --root results/wetting
```

Use fresh result roots. The scripts preserve executable hashes, generated inputs
and run receipts, and reject overwritten cases. Long droplet equilibrations
need the continuation procedures documented in the
[study README](research/drlab-publication/README.md) and supplement.
Recorded measurements and figure generation are available with the
[paper sources](research/drlab-publication/paper/).

## Input And Restart Conventions

- Geometry is row-major unsigned bytes, x fastest. `.prm` dimensions are in
  reversed order: `z y x` for 3D, `y x` for 2D. Every voxel label must be defined.
- The historical key is `rexlationTime` (spelling retained). `latticeSpeed 1`
  is required. `GV` specifies force density along a non-rest lattice direction.
- SCSP reports reference-density velocity; ColorFluid uses local total density.
  These coincide in the incompressible limit. Account for density when comparing
  their pure-phase fields on a rock with a nonuniform pressure field.
- `interfaceThickness` specifies recoloring beta, not the interface width in
  nodes. Solid `SD` label values specify contact angle through phase 1; the
  historical global `contactAngle` value is inactive.
- `wallGradient "legacy"` is the default; `"renormalized"` enables weighted
  moment reconstruction. `viscosityInterpolation "arithmetic"` is the default;
  `"harmonic"` selects harmonic kinematic viscosity.
- `wallGradient "tanh"` enables experimental local profile continuation;
  it does not uniformly improve the measured controls. `capillaryModel "csf"`
  selects the experimental curvature/Guo force instead of the default
  `"perturbation"`. CSF accepts only equal reference densities, standard
  pressure weights and `latticeSpeed 1`. Closed VD/GV flow uses SD/BB walls.
  Pressure CSF additionally requires the guarded flat D2Q9 ghost/flux
  combination described below. Velocity reservoirs and selective films
  remain unsupported. New quantitative capillary
  measurements use equal viscosity, tau=1; moving/curved 3D contact-line
  accuracy remains unverified.
- `wallRecoloring "tangential"` selects an experimental D2Q9 flat-wall
  correction with `wallGradient "renormalized"`. At 90 degrees it preserves
  the bulk discrete planar phase profile using a lattice-weight-derived
  factor and the true tangent-gradient mask. Equal reference densities and
  standard weights are required; corners, staircase walls and open boundaries
  are rejected. Nonneutral tests are exploratory and show larger pressure
  bias. This option is not validated for digital-rock wall geometry.
- `curvatureStencil "gradient"` is an optional CSF discretization of the same
  continuum force; `"normal"` remains the default. The gradient form lowers
  neutral residual speed but worsens the tested bulk Laplace pressure error.
- `phaseBoundary "flux"` selects a conservative relative-phase transfer at
  flat D2Q9 walls; `"bounce"` remains the default. It restores the complete
  neutral bulk phase recurrence and preserves pure-phase no-slip reflection.
  It requires `wallRecoloring "legacy"`, equal reference densities, standard
  weights and a supported reconstructed wall gradient. Open CSF requires
  the separately validated pure-collar pressure configuration below.
  Transport consistency alone does not improve nonneutral pressure accuracy.
- `wallGradient "tangent"` preserves the central wall-tangent derivative;
  `"halfway"` extrapolates the inward normal component to the physical wall;
  `"ghost"` uses the published flat-wall virtual-color formula of Liu et al.
  (2015), DOI 10.1103/PhysRevE.92.033306, without a second angle rotation.
  These experimental controls require complete flat D2Q9 stencils and
  nonsingular angles; halfway/ghost also require an interior fluid layer.
  The ghost option retains LBFlow's other operators, so it is not the full
  published model. Measured pressure improvement retains nonneutral flow
  and two-step oscillations. Curved digital-rock walls remain unsupported.
- `recoloringModel "projection"` uses the published link projection in
  place of the recovered `"cosine"` normalization. It requires equal-density,
  standard-weight D2Q9 and cannot combine with the cosine-derived tangent
  multiplier. Exact axial transport matching at cosine beta 0.7 gives
  projection beta 0.6316582489435276; equal beta changes the interface width.
- `mrtSpectrum "liu2015"` selects the published D2Q9 rates with matching
  Guo forcing; `"legacy"` remains the default. It requires color MRT and
  `capillaryModel "csf"`. Its fixed odd-moment rate changes the hydraulic
  wall offset; do not inherit the legacy pure-channel reference.
- Pressure CSF is experimental: use two complete opposite inward axial PB
  planes, `wallGradient "ghost"`, `phaseBoundary "flux"`, flat parallel
  walls and no external body force. The prescribed phase must remain pure
  within 1e-9 at PB and two interior strips; the first strip must have
  gradient below 1e-9. The run fails when an interface reaches this collar.
  `LBFlowBoundaryBudget` records interior component mass and exact kinetic
  inflow/outflow, with a new baseline on restart. Long physical validation
  covers theta60, equal tau1 and H16/H32. These are flat-slit drift controls,
  not validated curved-pore entry thresholds or contact-line laws.
- `wallSlipLength 1.0` prescribes a lattice-unit hydraulic Navier slip length;
  zero is the default. This research control currently requires neutral flat
  D2Q9 SD/BB walls, equal tau=1, standard equal-density color flow,
  projection recoloring, legacy MRT and CSF ghost/flux. Closed VD/GV labels
  are admitted; pressure reservoirs and nonneutral/curved walls are rejected.
  Forward specular populations and complementary component transfer preserve
  the neutral phase recurrence. Pure-channel calibration, rotation, MPI and
  restart are checked; molecular slip and rock contact-line laws are not
  inferred. Under grid refinement, convert physical slip to lattice units
  using `lambda_lattice=lambda_physical/h`.
- `reportDir 0` requests a prescribed run duration. Automatic convergence
  requires three stable velocity/color reports and applicable scalar diagnostics,
  including changes from the adjacent time step. Even reporting intervals
  can conceal period-two lattice modes; adjacent checks prevent false stopping.
  Legacy K is a reference-viscosity mobility diagnostic, not verified phase
  permeability. Axisymmetric `AB` boundaries and Shan--Chen are unsupported.
- Restart populations are post-collision. Keep the `.rst.meta` sidecar, geometry
  and physical configuration. Set `startTime` to the checkpoint time and provide
  `outPath.TIMESTEP.rst` and its sidecar. Old unlabeled checkpoints are rejected.
  Convergence history resets on load; native-double checkpoints are not a
  cross-platform interchange format.

At completion, `LBFlowTiming` provides a JSON record of update counts, fluid
nodes, numerical time (including MPI/finite checks), field/diagnostic time,
checkpoint time and initialization. Durations are maxima across ranks; their
sum need not equal the maximum run wall time. Construction and launcher costs
are measured separately by the Python runner. Dedicated Linux timing scripts
pin one logical CPU per MPI rank, discard warm-ups and repeat fixed workloads.

## Author And License

Yu Zhou, Independent Researcher, Canada.
LBFlow code is licensed under the [MIT License](LICENSE). The Bentheimer data
retain CC BY 4.0, and the archived cJSON source retains its upstream copyright
and license notices. Use [CITATION.cff](CITATION.cff)
and record your exact simulation commit when citing the software.
