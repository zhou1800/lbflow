# LBFlow Solver

Configure and build from the repository root; see `../README.md`.
`inc/`, `src/` and `apps/LBSim/` contain the C++ solver and its numerical operators. `tests/numerical_audit.cpp` emits production operator
coefficients for independent checks. The additional test sources, image
utilities and helper scripts are preserved but are not part of the supported
build or quantitative validation.

Operator verification and source checksums are in
`../research/drlab-publication/source-manifest.json` and
`../research/drlab-publication/notes/numerical-changes.md`. Generated configuration
headers live in the build tree.
