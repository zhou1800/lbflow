# Small Simulation Examples

Generate either a 3D force-driven channel (8 x 10 x 4 voxels) or a 2D
two-phase drop (24 x 24 nodes). Each command creates a new directory containing
an unsigned-byte geometry and a DRLab parameter file. Existing directories
are rejected to protect previous results.

From the repository root, after building:

```sh
python3 examples/generate.py channel --output results/channel
(cd results/channel && mpiexec -n 2 ../../build/bin/LBFlow input.prm)
python3 examples/generate.py drop --output results/drop
(cd results/drop && ../../build/bin/LBFlow input.prm)
```

The solver resolves paths relative to its working directory. Both inputs use
relative paths, so run from the generated directory. Outputs are `field.*.vti`,
population checkpoints `field.*.rst` with `.rst.meta` sidecars, and text logs.
ParaView can display the VTI fields.

The channel uses halfway bounce-back walls at y=0.5 and y=8.5, periodic x/z,
unit density, kinematic viscosity 1/6 and force density 1e-6 in +x. Its steady
velocity is `u(y) = F (y - 0.5) (8.5 - y) / (2 nu)` at fluid node centers.
The small drop demonstrates input and output; quantitative capillary
verification uses the larger cases in `scripts/benchmarks.py`.

The old C prototype's JSON format is archived with that prototype. The
maintained solver reads the line-oriented `.prm` format shown in these examples.
