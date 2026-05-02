# dos-forward-models

> [!CAUTION]
> **WORK IN PROGRESS**: This repository is currently being reorganized for public release. Documentation and examples are subject to change.

Forward modeling of light propagation and sensitivity mapping for Diffuse Optical Spectroscopy (DOS) and Near-Infrared Spectroscopy (NIRS).

## Overview

This repository contains a comprehensive suite of MATLAB scripts for modeling how light travels through biological tissue and calculating the sensitivity (Jacobian) of various measurement types to changes in optical properties. It includes analytical models, sensitivity mapping tools, and adjoint Monte Carlo integrations.

## Repository Structure

### Source Code (`src/`)
- **`analytical/`**: Fundamental solutions to the Radiative Transport Equation and Diffusion Approximation for various geometries (homogeneous, spherical, two-layer).
- **`jacobian/`**: Tools for generating sensitivity maps (Jacobians) for complex (FD), continuous (CW), and temporal (TD) measurement domains.
- **`monte-carlo/`**: Adjoint Monte Carlo models for sensitivity mapping, integrated with MCXLAB.

### Examples (`examples/`)
- **`fluorescence_flow/`**: Simulation suite for dual-ratio fluorescence measurements in flow, including SNR modeling.
- **`pulse_oximetry/`**: Modeling the impact of LED spectral linewidth on pulse oximetry accuracy.

### Shared Data (`data/`)
- Consolidated datasets and Monte Carlo output files (e.g., `MCout_avg.mat`).

## Citations

If you use this toolkit in your research, please cite the relevant publications:

1.  **Dual-Slope Foundations:** Blaney, G., Sassaroli, A., Pham, T., Fernandez, C., & Fantini, S. (2019). Phase dual-slopes in frequency-domain near-infrared spectroscopy for enhanced sensitivity to brain tissue: First applications to human subjects. *Journal of Biophotonics*, 12(11), e201960018. [https://doi.org/10.1002/jbio.201960018](https://doi.org/10.1002/jbio.201960018)
2.  **Array Design & Theory:** Blaney, G., Sassaroli, A., & Fantini, S. (2020). Design of a source-detector array for dual-slope diffuse optical imaging. *Review of Scientific Instruments*, 91(11), 114102. [https://doi.org/10.1063/5.0015512](https://doi.org/10.1063/5.0015512)
3.  **Broadband Spectroscopy:** Blaney, G., Curtsmith, P., Sassaroli, A., Fernandez, C., & Fantini, S. (2021). Broadband absorption spectroscopy of heterogeneous biological tissue. *Applied Optics*, 60(25), 7552-7562. [https://doi.org/10.1364/AO.431013](https://doi.org/10.1364/AO.431013)

## Author
Developed by Giles Blaney, Ph.D.

---
*This repository is a reorganized and documented version of a personal codebase, performed by Gemini CLI.*
