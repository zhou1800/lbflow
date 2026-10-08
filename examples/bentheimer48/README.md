# Bentheimer Digital-Rock Example

These 48 × 48 × 48 crops are derived from the public, segmented 4.54 µm
Bentheimer image by Shah, Crawshaw and Boek,
[DOI 10.6084/m9.figshare.1200096.v1](https://doi.org/10.6084/m9.figshare.1200096.v1).
The data and their derivatives retain **CC BY 4.0**. See `DATA-LICENSE.txt`;
the repository's MIT license applies to its code. Data courtesy of Qatar
Carbonates and Carbon Storage Research Centre, Imperial College London.

`manifest.json` records source and crop checksums, axis order, crop coordinates,
voxel size, connected components and preprocessing. No isolated pores were
removed. Opposite faces are periodically connected by the simulator. This
small repeated specimen is a numerical example, not a representative volume
or a laboratory reference for the sandstone.

In `single-phase.raw`, 75 is pore and 255 is grain. In `two-phase.raw`, 75 is
initial phase 0 (periodic distance to the grain mask ≥ 3 original voxels),
76 is phase 1, and 255 is grain. The distance is an initialization aid,
measured between voxel centers; it is not a calibrated pore/throat radius.

From the repository root, after building and installing analysis dependencies:

```sh
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset examples/bentheimer48 --root results/rock-single --mode single --ranks 2
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset examples/bentheimer48 --root results/rock-two --mode two --ranks 2
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset examples/bentheimer48 --root results/rock-contracts --mode contracts
```

Each command needs a fresh result directory. The single-phase study tests
force linearity and both pure-phase viscosity reductions. The two-phase
study records evolving phase geometry and advective flux, with a zero-force
control. The default interface is about five lattice nodes wide, comparable
to many pores here. Curved-wall wetting accuracy and a steady relative-
permeability curve are **not established** by this example.

Recreate the distributed crops, or subdivide each original voxel twice along
each axis while preserving its phase label and binary box geometry:

```sh
python3 scripts/prepare_bentheimer.py --output results/bentheimer48 --size 48
python3 scripts/prepare_bentheimer.py --output results/bentheimer48-r2 --size 48 --refinement 2
python3 scripts/rock_workflow.py --solver build/bin/LBFlow --dataset results/bentheimer48-r2 --root results/rock-refinement --mode refinement --ranks 4
```

Source download and extraction are checksum checked. Subdivision adds no new
micro-CT information. The refinement workflow scales force, surface tension
and duration diffusively; it measures sensitivity without asserting an
asymptotic convergence regime.
