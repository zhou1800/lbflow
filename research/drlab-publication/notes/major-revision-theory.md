# Major Revision: A Predictive Numerical Resistance Mechanism

This note derives and interprets the discrete interface-kinetics and
capillary-resistance experiments. Numerical diagnostics do not establish
a new physical contact-line law.

## Selective Intervention

The isolated-interface expansion gives a numerical shape coefficient
`K_num = D h^2/4`. At the tested `D=1/6` this is `h^2/24` in the common base
units. Adding `-a K_num L_y^2(phi)` to the independent FV reference preserves
phase mass and its planar equilibrium and does not change leading diffusion,
hydrodynamics, prescribed slip or width. The fixed strengths are 0, 0.5 and 1;
none is fitted to a moving-front result. This transverse diagnostic is not
a new physical phase equation or the exact complete LBM phase stencil.

The wall/bulk partition uses the same outer neutral Laplacian with an
intermediate nonnegative mobility weight. Its two pieces sum to the uniform
operator, preserve mass and have nonpositive quadratic work. The prescribed
wall window has length two in the H16 units, one eighth of slit height;
it is not a measured microscopic contact-line length.

## Reduced Response

At low force write a paired interface perturbation as
`delta phi(x,y) = -p'(x) eta(y)`, retaining the independently specified planar
normal profile. Project conservative phase advection onto the integral of
the ascending interface. Linearizing capillary forcing and solving the
verified Stokes problem gives

`eta_t = g F + (sigma C_h - a K_num L_y^2) eta`.

Here `g` is the independently calculated pure-body-force response, and
`C_h` includes the given profile, width, hydrodynamic slip and reference
curvature/face sampling. No evolving front field, force audit or front
velocity is used to construct either quantity. For steady translation,

`(sigma C_h - a K_num L_y^2) eta + g F = U 1`.

Fix the displacement gauge with `sum(eta)=0`. Equivalently the left null
vector `ell_a` gives the phase-displacement mobility

`U/F = (ell_a^T g)/(ell_a^T 1)`.

The physical mean fluid velocity is also calculated from the resulting
capillary/body force. It remains a separate recorded observable. The small
normal-profile reduction error is measured against full nonlinear reference
trajectories, not assumed to vanish. A centered finite-difference repeat
checks the response-matrix construction. Rigid translation has zero
capillary force; its small finite-difference symmetry residue is recorded
before imposing that analytic identity.

At `a=0`, multiplying the capillary matrix by sigma changes the displacement
needed to balance force but not the reduced mobility. At nonzero `a`, the
ratio `K_num/sigma` changes the left null vector and therefore the mobility.
This predicts a numerical surface-tension sensitivity in the creeping-flow
response without fitting a mobility or an effective slip. Full FV can retain
small profile/time effects beyond this held-profile reduction.

With capillary time scale `mu H/sigma`, the corresponding dimensionless
numerical strength is `Xi = mu K_num/(sigma H^3)`. Width/height, slip/height
and spatial discretization remain additional inputs. This is not a universal
error threshold. Symmetric parabolic driving excites even transverse modes;
the first asymmetric wave alone does not characterize the resistance.
The wall Stokes operator also mixes those modes. Initial projected capillary
rates are recorded separately from full coupled eigenvalues.

## Evidence And Scope

Three-grid slip1 interventions bring mean speeds within about 1.02%, 0.45%
and 0.28% of LBM, compared with the original FV gaps 17.12%, 15.18% and
13.72%. Slip4 and surface tensions .0025/.01 are held-out controls. The
intervention coefficient comes from the earlier derivation. The reduced
model independently predicts both untreated and treated FV response; the regenerated
mean-speed errors are at most 0.00875% in the baseline and 0.0564% across
all 22 available nonlinear-reference comparisons. The retained-normal-profile
volume-shape error reaches 4.244% in FV; mean-speed accuracy does not imply
an exact shape reduction.

The continuous-source manufactured calculation exercises nonzero surface
tension, analytic moving interfaces and nonuniform Navier flow. Phase and
velocity errors converge near second order at H16/H32/H64 for slips 1 and 4.
Its C2 logistic-to-pure-phase joins remove distance cusps in the bulk. Log
reconstruction remains regularized; small tail overshoots are retained and
reported, not clipped. Manufactured consistency is not an unforced physical
moving-contact-line solution.

On the same audited production phase fields, complete actual-force Stokes
velocity differences decrease to about 0.093% at H48. The raw capillary-force
norm difference is about 2.08%, and the difference in the Stokes velocities
it drives is about 0.100%. These are distinct measurements: pressure-gradient
components and viscous filtering can make force and velocity norms behave
differently. The entire force difference is not assigned to a pressure gauge.

The primary contribution is a quantified numerical phase-kinetics/resistance
mechanism in a specified neutral equal-viscosity slit. Conservative Allen--Cahn
mapping, adjustable mobility, signed-distance profile preservation and Navier
mixing have prior work. Haghani et al. (2024) additionally treats density
contrast and hydrodynamic/interface decoupling; that model is not reproduced
here. The new results do not establish molecular wetting, a universal contact
line law, high density contrast or quantitative curved-rock displacement.
The small transient rock example remains a demonstration.
