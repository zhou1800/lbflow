# Numerical changes and evidence

Reference source revision: `1bef9bbe7047ca74280f5870f4787be4df12fe20`. All changes below are maintained in the repository-root `solver/` tree.

| Change | Contract | Independent evidence |
|---|---|---|
| Include `<cstdint>` and update CMake | Build with a current C++11 toolchain | GCC 13.3 / CMake 3.28.3 / Open MPI 4.1.6 build; no forced-header workaround |
| Local MRT relaxation; zero conserved rates | Pure second-phase viscosity follows its own parameter; capillary-coupled moments share local rate | 24 analytic channel cases; emitted matrix/moment checks; original versus repaired unequal-viscosity Laplace cases |
| Guo source with half-force physical velocity | Correct momentum increment and stress in SRT/MRT and D2/D3 | 16 uniform acceleration cases; production-source moment identities; analytic wall channels |
| Complete initial collision from incoming rest state | Integer time and outgoing populations use the same force convention | Uniform velocity at time 0/50/100; uninterrupted/restarted fields |
| Canonical checkpoint stage and sidecar | Every checkpoint starts with the same next operation; no stale append or unlabeled stage | Six field comparisons, including one-to-two-rank restart; length and metadata checks |
| Field-based stopping | Constant saturation and velocity cannot stop a changing color field | Translating-indicator regression; historical premature-stop counterexample and repaired pressure case |
| Finite population/density checks and MPI failure status | Invalid state cannot become convergence or a valid checkpoint | Invalid input and corrupt-checkpoint executions; nonzero failure codes |
| Optional harmonic viscosity | Separate finite-interface closure choice from collision correctness | Twelve layered-flow executions; independent sharp and measured-diffuse references |
| Optional wall moment reconstruction | Recover locally linear gradients from the actual available stencil | Algebraic linear-field check; 30 wall configurations sharing one bulk kernel |

The uniform-viscosity MRT fast path caches its collision and force matrices. Unequal viscosities use local moment-space rates.
Perturbation and recoloring remain the default recovered mechanism; the capillary perturbation is not multiplied by the external Guo source prefactor.
The wall projection and wall normal remain common in the gradient ablation. Angle improvement is assessed alongside spurious velocity and fitting sensitivity.

The original native-double population layout is retained with explicit outgoing-state metadata. Restart histories reset conservatively; identical automatic stopping time across a restart is not claimed.
The sidecar does not identify every physical parameter. Run receipts contain the intended configuration and geometry hashes.
The executable workflow uses one configuration per MPI process lifetime; reuse of its historical function-static caches across multiple simulator objects is unverified.

Completed quantitative coverage is equal density, standard pressure weights, low speed and relaxation times 0.6--1.4. Wider density ratios, arbitrary EOS, curved/3D wetting, hysteresis, film boundaries and dynamic contact lines require additional work before new claims.

## Optional Capillary Research Operators

The continued wall study derives an optional `wallRecoloring "tangential"`
from the discrete D2Q9 segregation/diffusion ratio. It requires the
renormalized gradient and preserves the bulk neutral-plane fixed point.
At a neutral Cartesian wall it also removes the false raw normal derivative
before projection, so degenerate extrema retain the bulk active mask.
Corners, staircase walls, open reservoirs and nonstandard reference weights
are rejected. Nonneutral tests show larger pressure bias; this is not a
validated general wetting treatment.

The optional `curvatureStencil "gradient"` differentiates the unnormalized
color gradient. Its continuous force agrees with the default normal form,
but its matched bulk-drop pressure errors are larger. Both CSF forms require
unit lattice speed for the implemented Guo source. Scientific results use
the frozen `f3b646a...` kernel; reconstruction patches and exact source hashes
are in `run-kernels/`. The final unsupported-speed guard does not change
valid unit-speed fields, as three bitwise comparisons verify.

Convergence now checks adjacent physical velocity and phase values as well
as report-to-report changes. Even report intervals can alias the established
period-two checkerboard mode. Full nonlinear equilibrium seeds with SCSP and
color, SRT and MRT, exercise the production stopping path; stationary controls
still stop, while oscillatory controls complete their prescribed duration.
See `wall-study-validation.json` and `paper/data/wall-consistency.json`.

`wallGradient "tanh"` fits locally affine `atanh(chi)` on the available fluid
stencil, continues it over a virtual full stencil, and applies the common
geometric rotation. Color is clipped to `[-1+1e-8,1-1e-8]` for this transform.
This local reconstruction is inspired by equilibrium-profile continuation;
it is not an exact implementation of Bhattacharya et al.'s global ghost-node
scheme. Its H16 pilot reduces some nonneutral residual velocities but
increases the neutral-wall velocity. It is experimental, not a new default.

`capillaryModel "csf"` implements the established curvature surface force
`sigma*kappa*grad(chi)/2` through the existing Guo SRT/MRT operators, with
`kappa=-(I-n*n):grad(n)`. It exchanges a normalized-interface field and
reconstructs its derivatives from fluid neighbors at walls. The same total
force enters physical half-force velocity and collision. Fresh incoming
populations subtract the initial capillary half-force; checkpoints retain
their outgoing stage. This follows Lishchuk/Guo as used by Sedahmed and Coelho
(2024); our retained wall/color operators differ from their complete scheme.

The option rejects unequal reference densities, nonstandard pressure weights,
open reservoirs and selective films. Direct slit and bulk-drop measurements
use equal viscosity, tau=1. CTest additionally verifies SRT/MRT, D2Q9/D3Q19,
normal MPI exchange, initial bulk-plane equilibrium, uniform acceleration and
bitwise restart. These contracts do not transfer high-density-ratio or curved
3D wetting validation to this option. The old wall modes/default perturbation
retain bitwise-identical fields in four comparisons to the frozen baseline.

The new measurements show lower residual velocities with CSF and continuing
finite-width pressure bias. Zero-sigma controls retain a wall-color profile
distortion while their velocity is at roundoff. Neither a changed gradient
nor a force substitution proves consistent wall segregation or a new model.

## CPC Application And Pressure Reporting

The CPC extension adds independent pore-flux and pressure-entry controls,
public Bentheimer preparation, explicit phase advective flux, closed-mass and
late-field audits, and workstation profiling. Original archive files remain
unchanged. A density-normalized pure-phase comparison preserves SCSP's
reference-density velocity convention instead of treating it as an error.

SCSP now displays the standard pressure difference, Delta rho / 3, in its
pressure log. Previously the factor three was folded into the K prefactor
while the density difference was labelled pressure. The K formula is
algebraically unchanged at unit reference density. The pressure-driven slit
regression gives 19.00215 micron^2 against the independently sampled analytic
19.0 micron^2 (relative discrepancy 1.13e-4); its pressure discrepancy is
9.85e-5. No collision, force or wetting operator was changed by this log fix.

Profiling reports numerical work, initialization, fields/diagnostics and
checkpoint I/O separately, along with completed updates and physical fluid
voxel counts. Numerical work includes communication and finite checks.
Category maxima across ranks are not added to infer a wall duration.
