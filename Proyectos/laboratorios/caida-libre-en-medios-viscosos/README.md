# Caída libre en medios viscosos

Experimental physics project focused on the fall of spheres through a viscous fluid and the estimation of terminal velocity and viscosity using Stokes' law.

The **final academic deliverable for this experiment is a poster**, not a paper. The repository therefore keeps the poster inside this project together with the experimental data, Tracker files, videos, analysis code, results and figures.

## Experimental workflow

1. Record the fall of spheres of different sizes through the viscous medium.
2. Track the motion using Tracker.
3. Extract and organize position/velocity data.
4. Estimate terminal velocity for each sphere.
5. Use Stokes' law to estimate viscosity and check the flow regime with the Reynolds number.
6. Compare the behavior of the different spheres and prepare the final poster.

## Contents

- [`docs/poster.pdf`](./docs/poster.pdf) — final academic poster.
- [`data/experimental-data.xlsx`](./data/experimental-data.xlsx) — experimental workbook.
- [`data/processed/`](./data/processed/) — processed terminal-velocity and sphere-summary results.
- [`tracker/`](./tracker/) — Tracker project files for the three recorded spheres.
- [`videos/`](./videos/) — original experimental recordings used by Tracker.
- [`analysis/`](./analysis/) — Python scripts and saved analysis output.
- [`figures/`](./figures/) — figures generated from the analysis.

## Analysis notes

The available scripts implement velocity fitting and Stokes-law calculations. The saved analysis output includes terminal-velocity estimates, viscosity calculations and Reynolds-number checks for three sphere sizes. These files are preserved as part of the reproducible workflow; the poster is the final academic presentation of the experiment.

## Authors

**Paul Berdichewsky** · **Rafael Williams**  
Universidad Técnica Federico Santa María (UTFSM)
