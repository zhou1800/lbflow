# Numerical changes and evidence

Baseline: archived DRLab 0.2.2.1, revision `1bef9bbe7047ca74280f5870f4787be4df12fe20`. All changes below are maintained in the repository-root `solver/` tree.

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
Perturbation and recoloring remain the recovered color-gradient mechanism; the capillary perturbation is not multiplied by the external Guo source prefactor.
The wall projection and wall normal remain common in the gradient ablation. Angle improvement is assessed alongside spurious velocity and fitting sensitivity.

The original native-double population layout is retained with explicit outgoing-state metadata. Restart histories reset conservatively; identical automatic stopping time across a restart is not claimed.
The sidecar does not identify every physical parameter. Run receipts contain the intended configuration and geometry hashes.
The executable workflow uses one configuration per MPI process lifetime; reuse of its historical function-static caches across multiple simulator objects is unverified.

Completed quantitative coverage is equal density, standard pressure weights, low speed and relaxation times 0.6--1.4. Wider density ratios, arbitrary EOS, curved/3D wetting, hysteresis, film boundaries and dynamic contact lines require additional work before new claims.
