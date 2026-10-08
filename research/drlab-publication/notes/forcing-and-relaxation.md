# Forcing And Relaxation Contract

This repair uses lattice units `dx=dt=c=1` and the existing standard D2Q9/D3Q19 velocity weights.
The bulk pressure weights may depend on alpha; the velocity-polynomial and force quadratures remain the standard weights.
Claims outside the validated density/alpha range require additional equivalent-equation analysis.

For the local force density F, define physical velocity by `rho*u = j + F/2` (reference rho for the existing incompressible SCSP branch).
The Guo quadrature is `G_i = w_i [3(e_i-u) + 9(e_i.u)e_i].F`.
Its moments are `sum G_i=0`, `sum e_i G_i=F`, and `sum e_i e_i G_i=uF+Fu`.
The SRT source is `(1-1/(2*tau))*G`.
The MRT source is `M^-1 (I-S/2) M G`, with local `S=S(tau(x))` for the color model.
Equilibrium momentum differs from the pre-collision momentum by F/2.
Consequently equilibrium relaxation contributes `S_j*F/2` and the source contributes `(1-S_j/2)*F`; their sum is exactly F.
This also holds for the historical D2Q9 conserved-mode rates of one and D3Q19 rates of zero.

The executable streams first, measures macroscopic fields, then collides.
A checkpoint therefore stores a post-collision distribution.
At time zero with requested physical u=0, initialize the post-collision momentum to F/2 using `f_i <- f_i + G_i(u=0)/2`.
After the first stream, physical velocity is F/rho; uniform periodic acceleration must give `u(t)=t*F/rho`.
The time-zero output remains the prescribed physical zero velocity.
Every restart, including convergence-triggered restart, is saved after collision.

Keep the original no-force equilibrium, perturbation and recoloring formulas.
The color perturbation uses the same local inverse relaxation time as its stress relaxation.
Validate surface tension numerically rather than asserting that changing MRT rates alone proves capillary consistency.
