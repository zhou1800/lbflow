# Publication Study Contract

This study develops a defensible manuscript from the author's DRLab solver and geometric wetting implementation.
The historical working copies remain unchanged.
The imported solver, numerical experiments, equations, source comparisons and manuscript are separate reviewable artifacts.

## Source And Scope

- Baseline: recovered the reference solver, HEAD `1bef9bbe7047ca74280f5870f4787be4df12fe20`; file hashes are in `source-manifest.json`.
- Source of truth for initial results: the preserved `drlab-assessment` artifacts from 2026-10-07.
- Maintained implementation: repository-root `solver/`; numerical changes are documented individually.
- Published measurement summaries and figures: `paper/data/` and `paper/figures/`. Raw outputs are generated in user-selected result directories and kept outside version control.
- Main application: pore-scale immiscible two-phase flow in porous media and micro-CT digital rock. The numerical scope remains bounded by completed validation.

## Numerical Outcomes

1. MRT relaxation follows the local two-phase relaxation time, including the pure second-phase limit.
2. Force integration produces the prescribed net momentum increment in both SRT and MRT and both supported lattices.
3. Reported velocity uses the same force convention as the collision operator.
4. Periodic constant-force acceleration and pressure/force-driven channels agree with their analytic references under refinement.
5. Closed-domain phase masses are conserved; phase geometry and parasitic velocity are diagnosed separately.
6. Convergence requires finite numerical state and velocity convergence as well as the requested report diagnostics.
7. Restart files always represent the same post-collision state, including a convergence-triggered restart.
8. Invalid numerical state and exceptions terminate with a failure status; no numerical failure is reported as convergence.

## Wetting Comparison

The historical truncated-fluid stencil is the control.
Any proposed gradient reconstruction or projection changes only the wetting treatment on a common verified bulk solver.
Compare against directly relevant published schemes and current preprints before claiming novelty.
Record the physical bounce-back wall location, phase through which the angle is measured, measured interface thickness, fitting exclusion height, initialization and convergence criterion.
Use at least three resolutions for a convergence claim.
Control nondimensional parameters when changing domain resolution and distinguish spatial resolution from diffuse-interface thickness.
Record negative results as well as improvements.

## Acceptance Evidence

- Pure-phase viscosity response at two unequal relaxation times with analytical channel checks.
- Uniform periodic acceleration: mean momentum/physical velocity versus prescribed force and time.
- Planar interface and curved-drop phase mass, pressure jump and parasitic velocity.
- Wall-drop angle versus target at several angles and resolutions; fitting sensitivity and late-time drift.
- MPI global-field agreement and uninterrupted versus restarted evolution.
- Source-level change log, run receipts, executable hashes and reproducible inputs.

## Build And Validation Contract

The repository-root CMake project builds the maintained C++11/MPI solver as `build/bin/LBFlow`.
The incomplete C rewrite and its old tooling remain under `legacy/c-prototype/`.
There is one maintained copy of the solver and of the experiment scripts.
The supported build works out of tree without changing source files or downloading dependencies during configuration.
The installation contains a working simulator with a documented generated small-volume example.
CTest checks actual operator moments, physical forcing, two-rank agreement, restart continuity and invalid-state exits.
GitHub Actions repeats these checks using system MPI and Python dependencies.

## CPC And Porous-Media Extension

The main paper explains the research use, algorithms, software structure and independently measured operating range, with a CPC program summary. Detailed historical repairs remain traceable but do not substitute for an application.

1. Define Darcy velocity, phase-weighted advective flux, saturation and reference-viscosity conventions explicitly. Verify phase mobility against analytical isolated parallel pore channels at unequal phase viscosities and several resolutions before reporting it on a rock.
2. Add a capillary pore-entry or displacement benchmark with stated wall location, contact angle convention, pressure/forcing, measured interface width and an analytical reference where available. Record failures and finite-interface sensitivity rather than fitting an assumed result.
3. Use a publicly licensed, DOI-identified segmented digital rock. Preserve original metadata, file hashes, voxel size, crop coordinates, orientation, preprocessing and pore connectivity. Separate a finite cropped demonstration from a representative-volume or laboratory validation claim.
4. Demonstrate single-phase permeability and two-phase pore redistribution or displacement, finite-state/mass checks, late-time changes and MPI/restart agreement on this geometry. Report relative permeability only to the extent independently validated; do not reuse legacy K as a phase permeability.
5. Measure computation and I/O separately, with fixed inputs, repetitions, hardware/affinity/compiler/rank details and identical fields where relevant. A workstation timing experiment is not a cluster or GPU scaling result.
6. Compare scope and conventions with directly relevant pore-scale software and publications, especially LBPM and current OpenLB. Do not claim a runtime comparison unless both implementations were actually measured with a defensible matched configuration.
7. Rewrite the manuscript and supplement around the selected application; update references, figures, program metadata and packages from completed evidence. Funding, competing interests and affiliation country must come from the author.

## Physics Research Extension

The research question is how finite interface thickness and wall treatment affect
viscous transport and capillary trapping in resolved pore geometries.
The existing program manuscript remains a measured baseline; a new physical
claim requires a mechanism, a controlled test and a comparison with recent work.

1. Derive the continuum hydraulic response of a variable-viscosity slit from
   momentum balance. Separate diffuse constitutive bias, interface-profile
   changes and residual lattice error. Test both viscosity closures at a 9:1
   ratio, three pore widths and three recoloring parameters. Calibrate the
   bulk interface separately. Preserve low-speed diffusive force and surface
   tension scaling. Do not call these coupled sharp-interface/grid sequences
   a fixed-interface grid convergence study.
2. Quantify the wall stencil's directional and magnitude errors in the
   capillary stress. Distinguish an exactly reconstructed linear gradient from
   a nonlinear contact-line or thermodynamic consistency result.
3. Treat capillary-threshold shifts, artificial pinning and residual trapping
   as hypotheses until controlled invasion and retraction sequences resolve
   their pressure intervals and rate dependence. Use a specified throat
   length for confinement, not a rock distance-transform median.
4. Compare recent viscosity, wetting, interface-materiality and porous-medium
   studies before making a novelty claim. Established closures, interface
   thickness sensitivity and familiar trapping mechanisms are not new by
   themselves. AI is optional and is not part of the present research scope.
5. Put implementation history and parallel checks in support of the physical
   evidence. Update the paper only with completed derivations and measured
   results. Keep a bounded next-stage plan for unmeasured pore-throat physics.

## Review-Driven Continuation

The author authorizes further physical research and model improvements when
needed. Scientific peer review uses domain evidence, not software-document
templates. The review of source `422b277` identifies missing directly relevant
diffuse-interface analyses, reference-reconstruction uncertainty, and an
unclosed connection to quantitative two-phase pore physics.

1. Compare with Lai et al. (2025), Fuster and Mimoh (2024), and Reder et al.
   (2024). Retain the continuum analysis as an explanatory reference and
   establish additional significance through measured physical outcomes.
2. Record linear and shape-preserving reconstruction sensitivity for the
   eighteen existing layered cases. Their residual is a discrepancy against
   a specified continuum reconstruction; do not call it isolated lattice
   error without an independent uncertainty bound.
3. Measure closed slit menisci with direct bulk pressure differences, signed
   Young-Laplace references, phase mass and late geometry changes. Include
   neutral wetting as a planar mechanical-equilibrium control. Resolve
   static pressure before inferring an open-boundary entry threshold.
4. Diagnose the wall representation before introducing an optional correction.
   A correction must preserve the established defaults, bulk force/capillary
   operators, source provenance and replay contracts. Compare pressure,
   geometry and parasitic flow, including negative results. Published tanh
   profile continuation is a baseline idea, not evidence of novelty.
5. Use a common calibrated surface tension, interface width and wall location.
   Establish finite-time and resolution sensitivity before new physical
   predictions or journal submission. No new theory, universal correction or
   submission-readiness claim follows merely from adding a model option.

### Experimental Nonlinear Wall Reconstruction

The static neutral-wall control shows persistent flow with both existing
wall gradients. An optional `wallGradient "tanh"` examines one specific
mechanism: linear reconstruction of a nonlinear equilibrium color profile.
It reconstructs the locally affine `atanh(chi)` field using the available
fluid-stencil moment tensor, continues the inferred profile over a virtual
full stencil, and applies the common geometrical angle projection. The
transform clips color to `[-1+1e-8,1-1e-8]`; this regularization is explicit.
Rank-deficient stencils remain unsupported. The legacy default and existing
renormalized operator retain their arithmetic and bulk implementation.

This is an experimental local reconstruction, not an implementation claim
for a published global ghost-node algorithm, a new bulk theory or a guarantee
of wall mechanical equilibrium. Test direct capillary pressure, neutral-wall
flow, mass, gradient-profile consistency and replay. Promote it only to the
scope supported by those comparisons; record any failure to improve flow.

### Experimental Curvature Surface Force

The nonlinear wall pilot does not remove the neutral-plane parasitic flow.
An optional `capillaryModel "csf"` replaces the perturbation increment with
the established continuum surface force `F_s = sigma*kappa*grad(chi)/2`,
where `kappa = -(I-n*n):grad(n)` and `n = grad(chi)/|grad(chi)|`.
This follows the Lishchuk/Guo force approach used in Sedahmed and Coelho
(2024), DOI 10.1063/5.0228835; it is not a new physical theory or an exact
replication of their complete wetting model. Normal derivatives use the
weighted fluid stencil, reconstructed at walls, with an additional MPI
normal-field exchange. Gradients below the existing recoloring cutoff have
zero normal and force. The geometric projection and chosen gradient remain
common to both capillary options.

The same combined external and capillary force enters Guo collision and
the physical half-force velocity. Fresh incoming populations include the
initial capillary half-force subtraction. Saved states remain outgoing.
The default `perturbation` arithmetic and archived source are preserved.
Initially admit only closed equal-density, standard-pressure-weight color
flow with VD/GV labels and latticeSpeed=1, as required by the existing Guo
source coefficients; reject unsupported speed units explicitly. Pressure/velocity reservoirs and selective films are
unsupported for this experimental option and must be rejected explicitly.

Verify zero-force plane equilibrium, signed slit pressure, Laplace pressure
of bulk drops, mass, both collision operators and lattices, MPI exchange,
restart continuity and unchanged default fields. Report pressure/speed and
late-time sensitivity together; lower spurious velocity alone does not
establish an accurate or universally superior wetting model.

### Discrete Wall Segregation Balance

The author requests continued resolution of the residual wall inconsistency.
First isolate the zero-force, unit-density phase dynamics with an independent
D2Q9 recoloring/streaming recurrence. Reproduce the recorded zero-sigma wall
profile before inferring a correction. Prior boundary-layer analyses, notably
Ginzburg, Roux and Silva (2015), DOI 10.1016/j.crme.2015.03.004, and the
zero-interfacial-force wetting approach of Yu et al. (2018), are direct
comparisons; do not claim the general boundary-layer observation as new.

For a neutral planar interface varying tangentially along a flat Cartesian
wall, write coordinate and diagonal equilibrium weights as a=1/9, b=1/36.
The bulk tangential diffusion coefficient is A=a+2b and the segregation
coefficient is B=a+sqrt(2)*b. The wall receives A_w=a+b and B_w=a+b/sqrt(2).
An optional `wallRecoloring "tangential"` multiplies recoloring only on the
two wall-parallel coordinate links by
`k=[(a+b)*(a+sqrt(2)*b)/(a+2b)-b/sqrt(2)]/a`.
The wall's segregation/diffusion ratio then matches the bulk. For an exactly
one-dimensional profile and a common active-gradient mask, its wall update
residual is `(A_w/A)` times the bulk residual. This preserves the bulk
discrete stationary profile, not the complete transient transport operator.

Opposite links share the same factor: the added component mass sums to zero,
and complementary recoloring leaves every total population unchanged.
The optional rule initially supports closed equal-reference-density,
standard-weight D2Q9 with flat Cartesian walls. Reject corners, curved or
staircase wall stencils, open reservoirs and selective films. Nonneutral
angles on supported flat walls are experimental tests of this extension;
the neutral-plane derivation is not proof of their accuracy. The default
recoloring and archived source remain unchanged. No pressure, velocity or
surface-tension fitting parameter is introduced.

An odd-length manufactured equilibrium exposes a further degenerate case:
at a profile extremum the tangential gradient can vanish while the truncated
nonlinear stencil creates a false normal derivative. The old rotation turns
this derivative into an active tangential recoloring direction. The optional
balanced rule therefore requires the renormalized wall stencil and, at a
neutral Cartesian wall, removes the raw normal component before computing
the projected magnitude and common cutoff. For a one-dimensional profile,
the reconstructed tangent then equals the bulk central derivative, including
its inactive mask. This is part of the neutral fixed-point correction, not
a changed cutoff or a general nonneutral gradient reconstruction. Retain the
failed odd-grid control in the research evidence and rerun it after the fix.

Verify independent phase-recurrence agreement, stationary-profile residuals,
zero/nonzero-sigma neutral controls, three beta values and slit widths,
surface-force versus perturbation interaction, component mass, MPI/restart,
invalid-stencil exits and unchanged defaults. Retain failed controls and
cutoff sensitivity. Only then advance to dynamic pressure-boundary/throat
tests; a static correction is not a contact-line mobility or trapping law.

### Curvature Discretization Control

The first correction eliminates the zero-sigma phase distortion, but normal
differentiation leaves residual flow after a step-profile initialization.
Compare the default `curvatureStencil "normal"` with an optional
`"gradient"` form of the same continuum force:
`F_s = -(sigma/2) n [(I-n*n):grad(g)]`, `g=grad(chi)`.
For a smooth unit normal, this follows the product identity
`(I-n*n):grad(g) = |g| (I-n*n):grad(n)`; it is not a new physical law.
The discrete forms are distinct and must be measured separately.
The gradient form differentiates neighboring unnormalized gradients, avoiding
division by weak neighboring gradients. It uses the same local normal,
wall moment reconstruction, force cutoff and Guo integration. Do not raise
the cutoff to hide currents or fit a surface-tension multiplier.

Test bulk Laplace pressure and neutral-wall phase/equilibrium controls before
nonneutral use; retain any loss of curved-interface accuracy. The default
normal form remains unchanged. New claims require mass, MPI/restart,
independent profile controls and a documented resolution/interface sequence.
The numerical-audit executable may inspect densities and capillary operators
from an identified outgoing checkpoint through `--capillary-state`. It does
not label post-collision momentum as physical velocity. Retain force and
gradient fields so weak-gradient and boundary defects can be distinguished.

### Adjacent-step Convergence and Staggered Modes

Force-state inspection and a lifted discrete equilibrium must distinguish
stationary force imbalance from an undamped temporal mode. An even reporting
interval can alias a period-two velocity: equal reported fields do not prove
stationarity. On a Cartesian lattice, the linearized equilibrium perturbation
`delta f_q=3*w_q*e_qx*A*(-1)^x` is unchanged by a momentum-conserving collision
and changes sign under streaming. Wall bounce-back also changes the sign of
this wall-parallel perturbation. This known checkerboard mechanism is not a
new mode or a new physical effect; compare Lallemand and Luo (2000) and the
two-phase analysis of Guo, Shi and Zheng (2011).

Save adjacent odd/even fields of the corrected slit, retain the raw speed,
measure the spatial Nyquist projection and adjacent-step difference, and
report an adjacent-step average only as a diagnostic. Do not filter the
solver fields or call this average a mechanical-equilibrium solution.
Manufactured phase equilibria must retain near-zero velocity under both
curvature forms, with their independently generated population seed retained.

Extend the existing convergence condition to require both report-to-report
and adjacent-step velocity stability, and (for color flow) both phase-field
conditions. Capture the previous physical velocity and phase only at a
report step, before computing that step's macroscopic state. Use the existing
absolute velocity floor and relative/phase tolerances; add no fitted damping.
Restart resets report history as before. Verify a full nonlinear equilibrium
checkerboard seed cannot falsely stop under even reports in SCSP and color,
SRT and MRT, while stationary controls still stop. This changes early-stop
decisions; fixed-duration collision, streaming and force arithmetic remain
unchanged.

### Conservative Wall Phase Flux

The next stage must distinguish a stationary-profile correction from a
transport boundary condition. Independently compare bounce-back with the
established nonlocal specular-forward scalar reflection (Ginzburg et al.,
2015). In a unit-density, zero-velocity neutral plane, the latter restores
the complete bulk tangential recurrence, not only its stationary ratio.

An optional `phaseBoundary "flux"` redistributes phase populations between
neighboring flat-wall fluid nodes. For a face joining wall nodes L and R,
let p/m be outgoing diagonal links with positive/negative tangent and normal
directed into the solid. Use
`J = f0_p(L)-f0_m(R) - [(phi(L)+phi(R))/2]*(fT_p(L)-fT_m(R))`.
Add J to the incoming diagonal at R and subtract J from its partner at L.
The second phase remains complementary to the original bounce-back total
population. Each face transfers zero net component mass; the total fluid
keeps its no-slip reflection. Pure phases produce J=0 even with nonuniform
total populations. At rest with constant total weights this is exactly
specular-forward scalar transport. The arithmetic face fraction is an
explicit numerical closure, not a derived general contact-line law.

Prepare outgoing phase/total population halos and their phase fraction
before streaming, including the rest population. This avoids relying on
the streaming buffer's stale rest direction or an uninitialized restart
macrostate. Initially support only closed equal-reference-density,
standard-weight, unit-speed D2Q9 on complete Cartesian flat-wall stencils.
Require `wallRecoloring "legacy"` so the flux correction does not also apply
the older multiplier. The neutral projection retains the true tangent
mask. Reject corners, staircase walls, selective films and reservoirs.

Also compare the established geometric tangent-based wall reconstruction
through `wallGradient "tangent"`: retain the central tangent derivative and
set the inward-normal derivative to `abs(g_t)*cot(theta)`. This differs from
rotating the full clipped magnitude. Flat D2Q9 stencils and nonsingular
angles are required; neither reconstruction is claimed as new. Related
color-conserving wetting closures of Hollis et al. (2011), Liu et al. (2015)
and Yu et al. (2018) require direct comparison before novelty claims.

Acceptance: independent face conservation and pure-phase tests; identical
one-dimensional neutral transient recurrence; full-solver zero/nonzero-sigma
neutral starts; signed nonneutral slit pressures and longer holds; actual
adjacent fields, component mass, both collision operators, rotation,
MPI/restart and unchanged defaults. Promote a general wetting claim only
if pressure, angle and velocity improve together across controlled widths.
Preserve failures. Dynamic reservoir/throat tests depend on these results.

### Contact-angle Location Control

The actual bounce-back wall is halfway between fluid and solid nodes, while
the current rotation imposes the desired normal at the first fluid center.
For a sharp circular meniscus this can replace H by H-1 and cause a leading
pressure shift of 1/(H-1); this is a testable geometric hypothesis, not a
surface-tension adjustment or an established explanation of every error.

Compare `wallGradient "halfway"` on flat D2Q9 walls. Extrapolate the unit
normal to the physical wall: `n_w=1.5*n_1-0.5*n_2` in its inward-normal
component. Set `c_1=(2*cos(theta)+c_2)/3`, use
`s_1=sqrt(1-c_1*c_1)` and preserve the first-row central tangent gradient
through `|g_1|=|g_t1|/s_1`. Obtain c_2 from the unprojected interior color
gradient. Where that interior gradient is below the established cutoff,
continue the target normal explicitly; reject singular reconstructed angles.
An interior fluid layer is required. Exchange the prepared bulk gradients
so a neighbor across an MPI boundary does not need a second halo layer.

This is a standard boundary-location/extrapolation control, related to the
halfway ghost treatments in existing geometric wetting work. It is not an
exact implementation of Liu et al. (2015). Check pressures, shape, actual
angle at the wall, neutral equilibrium, phase mass and MPI/restart against
the same physical spacing H. Do not fit an effective width to the results.

### Published Ghost-gradient Comparison

Implement `wallGradient "ghost"` as a direct comparison with the flat-wall
color continuation in Liu et al. (2015), equations 24--25. With this code's
phase/angle convention, set the virtual solid-node color to
`chi_ghost=chi_1-cot(theta)*abs(1.5*g_t1-0.5*g_t2)`.
Each row's tangent derivative uses a central difference. Exchange those
prepared derivatives and virtual colors at their associated fluid nodes,
then use all nine color-stencil entries, including diagonal ghost entries.
The gradient implicitly applies the angle at the physical halfway wall;
do not subsequently rotate it to impose the angle at the fluid center.
The virtual value is not a physical population and is not clipped or counted
as component mass. Require the same nonsingular flat-wall/interior-layer
conditions as the location control.

Compare local bounce-back and the conservative relative-phase flux using
unchanged recoloring weights, relaxation and curvature differentiation.
This is a faithful implementation of the cited scalar ghost formula under
the present remaining operators, not a replication of the entire 2015
model: that paper's link-length segregation weights and MRT spectrum differ.
Verify the ghost formula independently, signed pressure, neutral equilibrium,
angle/shape, speed, mass, MPI/restart and default preservation. Keep negative
results; do not infer a general curved-wall or dynamic contact-line result.

### Nonneutral Forcing And Published Core Operators

Separate transient lattice-mode excitation from maintained circulation.
For an even periodic Cartesian length with complete no-slip reflection,
let P be the spatially alternating physical momentum and G the corresponding
force projection. Streaming/reflection reverses the raw-momentum projection,
and local collision adds G. With physical half-force velocity, test the exact
relation `P(t+1)=-P(t)+[G(t+1)-G(t)]/2`. Independently read native outgoing
populations, physical VTI fields and the production force audit. Record
normalization, component, parity, storage stage and arithmetic residual.
Odd periodic lengths do not satisfy the even-periodic relation at the wrap.

Continue a held state with sigma set to zero as a diagnostic change of
physics. At equal densities, standard weights and equal tau, total-fluid
evolution then decouples from phase segregation. Compare decaying circulation
and surviving alternating momentum against unchanged-sigma continuations.
This is not a production velocity filter or a static-equilibrium result.
Compare even/odd domains and an exactly translated even geometry; retain
the small odd-domain geometry difference explicitly.

Provide independent optional controls `recoloringModel "projection"` and
`mrtSpectrum "liu2015"`, preserving defaults `"cosine"` and `"legacy"`.
Projection uses the published link-length-weighted segregation term
`beta*rho*phi*(1-phi)*w_i*(e_i dot n)` at standard weights. The D2Q9 MRT
spectrum in the existing row order is
`[1,1.63,1.14,1,1.92,1,1.92,1/tau,1/tau]`; all entries, including
conserved ones, act with the matching Guo factor and physical half-force
velocity. Admit only color, standard-weight equal-density unit-speed D2Q9;
the published spectrum requires MRT. Projection cannot combine with the
earlier cosine-derived tangential multiplier. Existing closed capillary
and flat-wall guards remain applicable.

At a neutral one-dimensional rest state, cosine segregation has
`B=a+sqrt(2)*b`, whereas projection has `B=A=a+2b`. An independent match
`beta_projection=beta_cosine*B/A` gives exactly the same bulk recurrence;
do not fit beta to pressure or velocity. Compare both the published beta
and this independently matched transport strength. Recalibrate the bulk
width for comparisons where exact one-dimensional matching does not apply.

Match every published operator that is explicitly specified, and identify
remaining boundary-normal/curvature continuation choices rather than
claiming a complete replica when details are unspecified. Verify force and
recoloring moments, pure-phase forcing, neutral transients, bulk Laplace,
signed nonneutral pressure/shape/speed, adjacent fields, translation,
MPI/restart and unchanged defaults. Existing checkerboard seeds must still
propagate their undamped mode: changing kinetic relaxation does not prove
that a conserved-momentum eigenmode has disappeared. Keep negative controls.
The fixed published MRT spectrum is supported with CSF only. Its energy
relaxation differs from the shear relaxation; do not reuse the legacy
perturbation-to-surface-tension mapping with that spectrum.
Pressure-reservoir and throat extensions depend on the measured gate.

For the published fixed odd-moment rate, also verify the hydraulic wall
response rather than inheriting the legacy magic-parameter result. At
low force in a pure D2Q9 channel, define
`Lambda=(tau-1/2)*(1/s_q-1/2)` and `C=4*Lambda/3-1/4`.
The midpoint Poiseuille profile is
`u=F*[y*(H-y)+C]/(2*nu)`, with y measured from the physical halfway wall.
Legacy rates give Lambda=3/16 and C=0; the published s_q=1.92 gives
Lambda=(tau-1/2)/48. Verify the constant offset and its decreasing relative
influence under refinement. A lower static current with a different MRT
spectrum is not proof of unchanged dynamic wall response.

### Guarded flat pressure-reservoir experiment

Extend CSF only for the explicit combination `wallGradient "ghost"`,
`phaseBoundary "flux"`, equal-density standard-weight D2Q9, flat Cartesian
walls and two opposite inward-facing axial PB planes. Reject velocity,
selective-film, interior pressure nodes, incomplete faces and externally
forced reservoir experiments. Each face prescribes a consistent phase and
positive density. Keep capillary force zero on PB nodes and require the
adjacent first fluid strip to be force-free; do not apply the uncorrected
Zou--He reconstruction through a forced interfacial layer.

At every incoming/outgoing state require the PB plane and two interior
strips to remain in their prescribed pure phase, within 1e-9 in phase
fraction. Also check the actual reconstructed color-gradient norm in the
first interior strip against the existing 1e-9 gradient cutoff. Fail if
the interface reaches this collar. Disable phase-only wall transfers on
every face touching a PB node, including the periodic wrap between opposite
PB planes. Preserve closed/default arithmetic.

Independently account for interior (VD/GV) component mass. Before each
stream, sum outgoing populations across PB-to-interior links as inflow and
interior-to-PB links as outflow. Exclude solid links, which return by
reflection. The interior mass change must equal the integrated net kinetic
flux; the boundary nodal velocity sum is not an exact replacement. Record
the initial interior mass, integrated flux, current mass and residual in
structured diagnostics. Reset this accounting baseline on restart and
state that reset; numerical populations themselves retain the exact stage.

Before reporting two-phase dynamics, check pure pressure-driven channels
against the selected MRT hydraulic response, two-direction mass budgets,
MPI/restart, and invalid geometry/collar guards. A closed held meniscus may
seed a single-front experiment by a recorded crop/translation and explicit
pure remote corridors while preserving every donor total population.
Surface tension and beta remain independently specified. Compare front
drift for pressure brackets, two durations and grid refinement. Call these
finite-time drift brackets, not a converged throat trapping or hysteresis
law. Adjacent-step averaging is diagnostic only; never filter solver state.
