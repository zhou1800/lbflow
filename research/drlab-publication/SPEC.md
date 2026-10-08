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
