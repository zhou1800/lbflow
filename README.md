# LBFlow

LBFlow is an MPI lattice Boltzmann simulator restored from Yu Zhou's Digital
Rock Lab (DRLab) code. The maintained C++11 implementation supports single-phase
flow and immiscible two-phase color-gradient flow on D2Q9 and D3Q19 lattices,
with SRT/MRT collision, voxel solids, geometric contact angles, VTI output and
population checkpoints.

The original C rewrite was an unfinished framework. It is preserved in
[`legacy/c-prototype/`](legacy/c-prototype/) with its tools and examples.
The runnable solver is now built from the repository root.

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

The restored solver corrects MRT phase viscosity and force integration,
physical velocity reporting, checkpoint time/stage semantics, field-based
convergence and numerical failure exits. The historical color-gradient
perturbation and recoloring remain the bulk mechanism. Optional harmonic
viscosity and wall-gradient reconstruction have separate measured tradeoffs.
See the [numerical change log](research/drlab-publication/notes/numerical-changes.md).

CTest exercises coefficients from the actual C++ implementation, 16 analytic
force cases, single/two-rank and restarted field agreement, and invalid parameter
and checkpoint failures. GitHub Actions builds with GCC and Clang and repeats
these checks and both small examples.

The extended study has 24 analytic channel cases, 12 layered-flow cases and 30
wall-drop configurations. Its completed quantitative scope is equal reference
densities, standard pressure weights, low speed and relaxation times 0.6--1.4.
Curved/3D wetting, moving contact lines and high density ratios need further
validation. MPI is implemented; CUDA and OpenMP acceleration are future work.

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
[paper sources](research/drlab-publication/paper/). The manuscript is a research
draft and has not been submitted or accepted.

## Input And Restart Conventions

- Geometry is row-major unsigned bytes, x fastest. `.prm` dimensions are in
  reversed order: `z y x` for 3D, `y x` for 2D. Every voxel label must be defined.
- The historical key is `rexlationTime` (spelling retained). `latticeSpeed 1`
  is required. `GV` specifies force density along a non-rest lattice direction.
- `interfaceThickness` specifies recoloring beta, not the interface width in
  nodes. Solid `SD` label values specify contact angle through phase 1; the
  historical global `contactAngle` value is inactive.
- `wallGradient "legacy"` is the default; `"renormalized"` enables weighted
  moment reconstruction. `viscosityInterpolation "arithmetic"` is the default;
  `"harmonic"` selects harmonic kinematic viscosity.
- `reportDir 0` requests a prescribed run duration. Automatic convergence
  requires three stable velocity/color reports and applicable scalar diagnostics.
  Legacy K is a reference-viscosity mobility diagnostic, not verified phase
  permeability. Axisymmetric `AB` boundaries and Shan--Chen are unsupported.
- Restart populations are post-collision. Keep the `.rst.meta` sidecar, geometry
  and physical configuration. Set `startTime` to the checkpoint time and provide
  `outPath.TIMESTEP.rst` and its sidecar. Old unlabeled checkpoints are rejected.
  Convergence history resets on load; native-double checkpoints are not a
  cross-platform interchange format.

## Author And License

Yu Zhou, Independent Researcher. Correspondence: zhouyu@live.com.
LBFlow is licensed under the [MIT License](LICENSE). The archived cJSON source
retains its upstream copyright and license notices. Use [CITATION.cff](CITATION.cff)
and record your exact simulation commit when citing the software.
