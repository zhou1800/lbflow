# Publication Study Contract

This study develops a defensible manuscript from the author's DRLab solver and geometric wetting implementation.
The historical working copies remain unchanged.
The imported solver, numerical experiments, equations, source comparisons and manuscript are separate reviewable artifacts.

## Source And Scope

- Baseline: recovered DRLab 0.2.2.1, HEAD `1bef9bbe7047ca74280f5870f4787be4df12fe20`; file hashes are in `source-manifest.json`.
- Source of truth for initial results: the preserved `drlab-assessment` artifacts from 2026-10-07.
- Maintained implementation: repository-root `solver/`; numerical changes are documented individually.
- Published measurement summaries and figures: `paper/data/` and `paper/figures/`. Raw outputs are generated in user-selected result directories and kept outside version control.
- First manuscript scope: immiscible two-phase color-gradient flow and geometrical wetting within the parameter ranges actually validated.

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

## Manuscript And External Tools

Prism is an authorized manuscript workspace.
Local LaTeX, references, figures and numerical data remain the reproducible source package.
ChatGPT review follows the explicitly invoked Tokimon ChatGPT skill, including verified project context, model/mode and paced observations.
Model output is advisory; equations, references and numerical claims require independent source or computational evidence.
No journal submission, author attestation, paid action or public release is performed before the concrete manuscript and submission materials are ready for the author's review.
The sole author is Yu Zhou, Independent Researcher, with correspondence at zhouyu@live.com.
The author has authorized migration to the existing GitHub repository and selected the MIT license.
Funding, competing interests and research adoption remain unresolved until the author supplies them.

## GitHub Migration Contract

The repository-root CMake project builds the maintained C++11/MPI solver as `build/bin/LBFlow`.
The incomplete C rewrite and its old tooling remain under `legacy/c-prototype/`.
There is one maintained copy of the solver and of the experiment scripts.
The supported build works out of tree without changing source files or downloading dependencies during configuration.
The installation contains a working simulator with a documented generated small-volume example.
CTest checks actual operator moments, physical forcing, two-rank agreement, restart continuity and invalid-state exits.
GitHub Actions repeats these checks using system MPI and Python dependencies.
The manuscript uses the confirmed author details; repository publication and journal submission remain separate actions.
