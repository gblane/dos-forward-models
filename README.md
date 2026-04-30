# dos-forward-models

Forward modeling of light propagation and sensitivity mapping for Diffuse Optical Spectroscopy (DOS) and Near-Infrared Spectroscopy (NIRS).

## Overview

This repository contains a comprehensive suite of MATLAB scripts for modeling how light travels through biological tissue and calculating the sensitivity (Jacobian) of various measurement types to changes in optical properties.

## Contents

### Sensitivity Mapping (`SenMaps/`)
Tools for generating sensitivity maps (Jacobians) across different measurement domains:
- **Complex (Frequency-Domain):** Models for Amplitude and Phase sensitivity (`complexFluence.m`, `complexReflectance.m`, `complexTotPathLen.m`).
- **Continuous (CW):** Sensitivity models for steady-state measurements (`continuousFluence.m`, `continuousReflectance.m`).
- **Temporal (Time-Domain):** Advanced models for time-resolved data, including moments (Mean Time, Variance) and gated measurements (`temporalReflectance.m`, `temporalKthMoment.m`, `temporalGatePartPathLen.m`).
- **Layered Models:** Sensitivity calculations specifically for two-layer media (`complexFluence2L.m`, `makeSenMaps2L.m`).

### Adjoint Monte Carlo (`SenMaps_MC/`)
Sensitivity calculations using Monte Carlo methods:
- `myMCXLAB_adjoint.m`: Integration with MCXLAB for adjoint sensitivity mapping.
- `temporalGatePathLen_MCadjoint.m`: Gated pathlength calculations using Monte Carlo.

### Analytical Forward Models (`abs_multiDist/` & `TwoLayer/`)
Fundamental solutions to the Radiative Transport Equation (under diffusion approximation):
- **Homogeneous & Spherical:** Solutions for standard geometries (`R_FD_forward.m`, `R_sphere_withPreCom.m`).
- **Two-Layer:** Efficient reflectance models for layered tissue (`R_2L_withPreCom.m`, `twoLayEffHomoOptProp.m`).

## Author
Developed by Giles Blaney, Ph.D.
