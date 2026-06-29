# Standardize Function Headers and Documentation

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Standardize MATLAB function headers and fix typos in the `dos-forward-models` repository.

**Architecture:** Use a systematic approach to read each `.m` file in `src/analytical/`, `src/jacobian/`, and `src/monte-carlo/`, extract existing documentation, and rewrite it using the specified template. Fix common typos during the process.

**Tech Stack:** MATLAB, Bash, Python (for scripting if needed), Git.

---

### Task 1: Fix common typos across all files

**Files:**
- All `.m` files in `src/analytical/`, `src/jacobian/`, `src/monte-carlo/`

- [ ] **Step 1: Run a global search and replace for common typos**

Run:
```bash
sed -i 's/distnace/distance/g' src/analytical/*.m src/jacobian/*.m src/monte-carlo/*.m
sed -i 's/frequecy/frequency/g' src/analytical/*.m src/jacobian/*.m src/monte-carlo/*.m
sed -i 's/corrdinates/coordinates/g' src/analytical/*.m src/jacobian/*.m src/monte-carlo/*.m
sed -i 's/abolsute/absolute/g' src/analytical/*.m src/jacobian/*.m src/monte-carlo/*.m
sed -i 's/measurment/measurement/g' src/analytical/*.m src/jacobian/*.m src/monte-carlo/*.m
```

- [ ] **Step 2: Commit typo fixes**

Run:
```bash
git add src/
git commit -m "docs: fix common typos in source files"
```

### Task 2: Standardize headers in `src/analytical/`

**Files:**
- `src/analytical/get_Rshpere_preCom.m`
- `src/analytical/PHI_TD_forward.m`
- `src/analytical/R_FD_forward.m`
- `src/analytical/R_FD_inf_forward.m`
- `src/analytical/R_sphere_withPreCom.m`
- `src/analytical/R_TD_forward.m`
- `src/analytical/get_R2L_preCom.m`
- `src/analytical/R_2L_withPreCom.m`
- `src/analytical/twoLayEffHomoOptProp.m`

- [ ] **Step 1: Update each file's header**
  For each file, extract:
  - Function signature
  - Brief description
  - Author and date
  - Inputs (with descriptions and units)
  - Outputs (with descriptions and units)
  
  Format into the new template:
  ```matlab
  function [outputs] = functionName(inputs)
  % functionName Brief one-sentence description of the function.
  %
  % [outputs] = functionName(inputs)
  %
  % Written by Giles Blaney, Ph.D. (Original date)
  %
  % Inputs:
  %   input1 - Description [units]
  %   ...
  %
  % Outputs:
  %   output1 - Description [units]
  ```
  Ensure units are in `[]` instead of `()`.

### Task 3: Standardize headers in `src/jacobian/`

**Files:**
- (List of 30 files from earlier search)

- [ ] **Step 1: Update each file's header using the template.**

### Task 4: Standardize headers in `src/monte-carlo/`

**Files:**
- `src/monte-carlo/continuousPathLen_MCadjoint.m`
- `src/monte-carlo/myMCXLAB_adjoint.m`
- `src/monte-carlo/temporalGatePathLen_MCadjoint.m`

- [ ] **Step 1: Update each file's header using the template.**

### Task 5: Final Review and Commit

- [ ] **Step 1: Verify all files have been updated.**
- [ ] **Step 2: Commit all changes.**

Run:
```bash
git add src/
git commit -m "docs: standardize function headers and improve source documentation"
```
