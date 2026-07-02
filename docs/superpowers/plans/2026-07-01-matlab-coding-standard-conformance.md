# MATLAB Coding-Standard Conformance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring all 57 `.m` files under `src/` and `examples/` into conformance with the §7 MATLAB coding standard without changing behavior, public names, signatures, or plotting output.

**Architecture:** A MATLAB regression harness captures golden outputs of every callable pure function *before* any edit. Files are then conformed subfolder-by-subfolder; each subfolder is gated by (a) output-equality regression against the golden baseline and (b) `checkcode` showing no new warning ids, then committed. Functions that cannot be exercised here (MCXLAB/GPU, complex struct inputs) fall back to `checkcode` + allowed-transformation review, reported explicitly.

**Tech Stack:** MATLAB (`/usr/local/bin/matlab`, `-batch` headless), `checkcode` (Code Analyzer), git.

## Global Constraints

Every task's work implicitly includes these, copied from the spec:

- **Behavior-preserving only.** No numeric output may change; regression gate uses `isequaln` (bit-identical).
- **No public renames.** Preserve every function file name and signature exactly.
- **No plotting rework.** Keep existing `figure/clf/plot/subplot` logic in examples.
- **Never alter** transpose operators (`'`, `.'`) or `%`-tokens inside string literals (`%d`, `%f` in `sprintf`/`fprintf`).
- **Formatting targets:** 4-space indent, no tabs, lines ≤ 120 chars, no trailing whitespace, no leading/trailing blank lines, single blank lines between logical blocks.
- **Strings:** single-quoted string *literals* → `"double-quoted"`; keep `'char'` only where a char vector is genuinely required.
- **Comments:** exactly one space after `%` (never touch `%%`, `%{`/`%}`, or `%` inside strings); simplify `%% ###…###` / `%% ---…---` banners to clean `%% Title`.
- **Operators:** spaces around `=`, relational, and logical operators; no spaces around `:`/`*`/`/`/`^`, inside brackets, in `Name=Value`, or after unary `+`/`-`/`~`.
- **Statements:** strict one statement per line (split `clear; home;`, `figure(101); clf;`, `subplot(...); hold on;`).
- **Approved skips:** keep cited/commented reference equations; do NOT convert `nargin`→`arguments`; do NOT sweep-rename variables (only touch a local var when a `checkcode` warning requires it).
- **checkcode gate:** no new warning ids per file vs. its pre-edit baseline.

---

### Task 1: Verification harness + golden baseline + checkcode-before

**Files:**
- Create: `verify/regressionSpecs.m`
- Create: `verify/captureBaseline.m`
- Create: `verify/verifyAgainstBaseline.m`
- Create: `verify/runCheckcode.m`
- Create: `verify/runBaseline.m`
- Create: `verify/runVerify.m`
- Baseline artifacts (git-ignored, scratchpad): `$SCRATCH/baseline.mat`, `$SCRATCH/checkcode_before.mat`

**Interfaces:**
- Produces: `regressionSpecs()` → N×4 cell `{name, fnHandle, argsCell, nOut}`; `captureBaseline(specs)` → struct array `(name, ok, out, err)`; `verifyAgainstBaseline(specs, baseline)` → struct array `(name, status)` where status ∈ `"PASS"|"FAIL"|"SKIP(baseline errored)"`; `runCheckcode(files)` → struct array `(file, ids, n)`.

- [ ] **Step 1: Confirm MATLAB launches headless**

Run:
```bash
matlab -batch "disp('MATLAB_OK'); disp(version)"
```
Expected: prints `MATLAB_OK` and a version string, exit code 0. If it fails (no license/headless), STOP and report; the pass downgrades to checkcode-only (record this and skip all regression steps).

- [ ] **Step 2: Write `verify/regressionSpecs.m`** (canonical inputs; conforms to the standard itself)

```matlab
function specs = regressionSpecs()
% regressionSpecs Canonical inputs for behavior-preservation regression tests.
%
% specs = regressionSpecs()
%
% Outputs:
%   specs - N-by-4 cell array {name, fnHandle, argsCell, nOut}. Functions that
%           cannot be constructed here are omitted and covered by checkcode only.
    op.nin = 1.4;
    op.nout = 1;
    op.musp = 1.2;
    op.mua = 0.01;

    rs = [0, 0, 0];
    rd = [25, 0, 0];
    r = [10, 0, 5];
    vol = 1;
    omega = 2*pi*140.625e6;
    z0 = 1/op.musp;
    c = 2.99792458e11;
    v = c/op.nin;
    tPs = linspace(10, 3000, 64);
    tSec = linspace(0.1e-9, 3e-9, 64);
    lam = (600:900).';

    specs = {
        "n2A",                  @n2A,                  {1.4, 1},                          1;
        "zeroOrdBesselRoots",   @zeroOrdBesselRoots,   {20},                              1;
        "LEDspec_func",         @LEDspec_func,         {(600:700).', 660, 30},            1;
        "tissueOptProps_func",  @tissueOptProps_func,  {lam},                             3;
        "Tslab",                @Tslab,                {25, 10, op, 100},                 1;
        "R_FD_forward",         @R_FD_forward,         {25, op.mua, op.musp, omega, v},   1;
        "R_FD_inf_forward",     @R_FD_inf_forward,     {25, op.mua, op.musp, omega, v},   1;
        "R_TD_forward",         @R_TD_forward,         {25, op.mua, op.musp, tSec, v, z0},1;
        "PHI_TD_forward",       @PHI_TD_forward,       {op.mua, op.musp, r, tSec, v},     1;
        "complexFluence",       @complexFluence,       {rs, rd, omega, op},               1;
        "complexReflectance",   @complexReflectance,   {rs, rd, omega, op},               1;
        "complexTotPathLen",    @complexTotPathLen,    {rs, rd, omega, op},               2;
        "complexPartPathLen",   @complexPartPathLen,   {rs, r, rd, vol, omega, op},       1;
        "continuousFluence",    @continuousFluence,    {rs, rd, op},                      1;
        "continuousReflectance",@continuousReflectance,{rs, rd, op},                      1;
        "continuousTotPathLen", @continuousTotPathLen, {rs, rd, op},                      2;
        "continuousPartPathLen",@continuousPartPathLen,{rs, r, rd, vol, op},              1;
        "temporalReflectance",  @temporalReflectance,  {rs, rd, tPs, op},                 1;
        "temporalFluence",      @temporalFluence,      {rs, rd, tPs, op},                 1;
        "temporalKthMoment",    @temporalKthMoment,    {rs, rd, 1, op},                   1;
        "temporalKthMomTotPathLen", @temporalKthMomTotPathLen, {rs, rd, 1, op},           1;
        "temporalVar",          @temporalVar,          {rs, rd, op},                      1;
        "temporalVarTotPathLen",@temporalVarTotPathLen,{rs, rd, op},                      1;
        "temporalGateTotPathLen",@temporalGateTotPathLen,{rs, rd, [500, 1500], op},       1;
        };
end
```

Note: functions absent here (`complexFluence2L`, `complexPartPathLen2L`, `complexTotPathLen2L`, `R_2L_withPreCom`, `get_R2L_preCom`, `get_Rshpere_preCom`, `R_sphere_withPreCom`, `twoLayEffHomoOptProp`, `complexFluence_Grad`, `complexReflectance_Grad`, `complexPartPathLen_Scat`, `complexTotPathLen_Scat`, `makeSenMaps*`, `makeS`, `sliceS`, `simMeasFromSenMap`, `temporalGatePartPathLen`, `temporalKthMomPartPathLen`, `temporalVarPartPathLen`, all `monte-carlo/*`) are **static-only** (checkcode + transformation review). Any spec row whose baseline call errors is auto-demoted to static and logged — a wrong input never masks a regression.

- [ ] **Step 3: Write `verify/captureBaseline.m`**

```matlab
function results = captureBaseline(specs)
% captureBaseline Evaluate each spec once, capturing outputs or the error.
%
% results = captureBaseline(specs)
    n = size(specs, 1);
    results = repmat(struct("name", "", "ok", false, "out", {{}}, "err", ""), n, 1);
    for i = 1:n
        name = specs{i, 1};
        fn = specs{i, 2};
        args = specs{i, 3};
        nOut = specs{i, 4};
        try
            out = cell(1, nOut);
            [out{:}] = fn(args{:});
            results(i) = struct("name", name, "ok", true, "out", {out}, "err", "");
        catch ME
            results(i) = struct("name", name, "ok", false, "out", {{}}, ...
                "err", string(ME.message));
        end
    end
end
```

- [ ] **Step 4: Write `verify/verifyAgainstBaseline.m`**

```matlab
function report = verifyAgainstBaseline(specs, baseline)
% verifyAgainstBaseline Recompute specs and compare to a golden baseline.
%
% report = verifyAgainstBaseline(specs, baseline)
    now = captureBaseline(specs);
    n = size(specs, 1);
    report = repmat(struct("name", "", "status", ""), n, 1);
    for i = 1:n
        name = specs{i, 1};
        if ~baseline(i).ok
            report(i) = struct("name", name, "status", "SKIP(baseline errored)");
        elseif ~now(i).ok
            report(i) = struct("name", name, "status", "FAIL(now errors)");
        elseif isequaln(now(i).out, baseline(i).out)
            report(i) = struct("name", name, "status", "PASS");
        else
            report(i) = struct("name", name, "status", "FAIL(output differs)");
        end
    end
end
```

- [ ] **Step 5: Write `verify/runCheckcode.m`**

```matlab
function summary = runCheckcode(files)
% runCheckcode Collect Code Analyzer message ids for each file.
%
% summary = runCheckcode(files)
    n = numel(files);
    summary = repmat(struct("file", "", "ids", strings(1, 0), "n", 0), n, 1);
    for i = 1:n
        msgs = checkcode(files(i), "-id");
        if isempty(msgs)
            ids = strings(1, 0);
        else
            ids = string({msgs.id});
        end
        summary(i) = struct("file", files(i), "ids", ids, "n", numel(msgs));
    end
end
```

- [ ] **Step 6: Write `verify/runBaseline.m`** (driver: adds paths per §2, captures baseline + checkcode-before, saves to scratchpad, restores path)

```matlab
function runBaseline(outDir)
% runBaseline Capture golden outputs and pre-edit checkcode summary.
%
% runBaseline(outDir)  -- outDir: folder to write baseline.mat / checkcode_before.mat
    srcPath = genpath("src");
    invPath = genpath(fullfile("..", "dos-inverse-models"));   % makeE
    myPath = genpath(fullfile("..", "my-matlab"));             % struct2pairs, saveFigure
    dsPath = genpath(fullfile("..", "dual-slope-toolkit"));
    addpath(srcPath, invPath, myPath, dsPath, fullfile(pwd, "verify"));
    cleanup = onCleanup(@() rmpath(srcPath, invPath, myPath, dsPath, ...
        fullfile(pwd, "verify")));

    specs = regressionSpecs();
    baseline = captureBaseline(specs);
    files = allSourceFiles();
    checkBefore = runCheckcode(files);

    save(fullfile(outDir, "baseline.mat"), "baseline", "specs");
    save(fullfile(outDir, "checkcode_before.mat"), "checkBefore", "files");

    okCount = sum([baseline.ok]);
    fprintf("Baseline captured: %d/%d functions regression-guarded.\n", ...
        okCount, numel(baseline));
    for i = 1:numel(baseline)
        if ~baseline(i).ok
            fprintf("  static-only: %-24s (%s)\n", baseline(i).name, baseline(i).err);
        end
    end
end

function files = allSourceFiles()
    d = [dir(fullfile("src", "**", "*.m")); dir(fullfile("examples", "**", "*.m"))];
    files = string(fullfile({d.folder}, {d.name})).';
end
```

- [ ] **Step 7: Write `verify/runVerify.m`** (driver: re-checks a file glob, errors nonzero on any regression/new-warning)

```matlab
function runVerify(outDir, fileGlob)
% runVerify Gate a set of edited files against the golden baseline.
%
% runVerify(outDir, fileGlob)  -- fileGlob e.g. "src/analytical" or "all"
    srcPath = genpath("src");
    invPath = genpath(fullfile("..", "dos-inverse-models"));
    myPath = genpath(fullfile("..", "my-matlab"));
    dsPath = genpath(fullfile("..", "dual-slope-toolkit"));
    addpath(srcPath, invPath, myPath, dsPath, fullfile(pwd, "verify"));
    cleanup = onCleanup(@() rmpath(srcPath, invPath, myPath, dsPath, ...
        fullfile(pwd, "verify")));

    S = load(fullfile(outDir, "baseline.mat"));
    B = load(fullfile(outDir, "checkcode_before.mat"));

    report = verifyAgainstBaseline(S.specs, S.baseline);
    nFail = 0;
    for i = 1:numel(report)
        fprintf("  %-26s %s\n", report(i).name, report(i).status);
        if startsWith(report(i).status, "FAIL")
            nFail = nFail + 1;
        end
    end

    if fileGlob == "all"
        sel = true(numel(B.files), 1);
    else
        sel = contains(B.files, fileGlob);
    end
    after = runCheckcode(B.files(sel));
    before = B.checkBefore(sel);
    nNew = 0;
    for i = 1:numel(after)
        newIds = setdiff(after(i).ids, before(i).ids);
        if ~isempty(newIds)
            nNew = nNew + 1;
            fprintf("  NEW WARNING %s: %s\n", after(i).file, strjoin(newIds, ", "));
        end
    end

    fprintf("Regression FAILs: %d | files with new warnings: %d\n", nFail, nNew);
    if nFail > 0 || nNew > 0
        error("runVerify:gate", "Verification gate failed (%d regressions, %d new-warning files).", ...
            nFail, nNew);
    end
    disp("VERIFY_OK");
end
```

- [ ] **Step 8: Run baseline capture on pristine code**

Run:
```bash
SCRATCH="/tmp/claude-1000/-home-giles-github-gblane-dos-forward-models/c7ff4041-2dc8-466f-99e6-4254cafe2342/scratchpad"
cd /home/giles/github/gblane/dos-forward-models
matlab -batch "runBaseline('$SCRATCH')"
```
Expected: `Baseline captured: K/24 functions regression-guarded.`, a list of any static-only demotions, exit 0. Record K and the demotions in the commit message.

- [ ] **Step 9: Commit the harness**

```bash
cd /home/giles/github/gblane/dos-forward-models
printf 'baseline.mat\ncheckcode_before.mat\n' >> .gitignore   # only if not already ignored elsewhere
git add verify/ .gitignore
git commit -m "Add verification harness for coding-standard conformance

Captures golden outputs (isequaln regression) + pre-edit checkcode ids so every
subsequent conformance edit is proven behavior-preserving. K/24 functions are
regression-guarded; the rest (2L/Scat/Grad/MC/makeSenMaps/makeS) are static-only.
Resume: re-run matlab -batch runBaseline before editing any pending subfolder.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Conform `src/analytical/` (14 files)

**Files (Modify):** `get_R2L_preCom.m`, `get_Rshpere_preCom.m`, `LEDspec_func.m`, `n2A.m`, `PHI_TD_forward.m`, `R_2L_withPreCom.m`, `R_FD_forward.m`, `R_FD_inf_forward.m`, `R_sphere_withPreCom.m`, `R_TD_forward.m`, `tissueOptProps_func.m`, `Tslab.m`, `twoLayEffHomoOptProp.m`, `zeroOrdBesselRoots.m`

**Interfaces:** Consumes the golden baseline from Task 1. Produces conformed files with identical behavior.

- [ ] **Step 1: Apply the Global-Constraints taxonomy to each file**

For every file, in order: strip trailing whitespace and blank-line padding; one space after `%` on comments (skip `%%`, `%{`/`%}`, and `%` inside strings); convert single-quoted string *literals* to double quotes (leave `'`/`.'` transposes untouched — e.g. `warning('...')`, `interp1(...,'linear','extrap')`, `error('...')`); add spaces around `=`/relational/logical operators and remove spaces around `:`/`*`/`/`/`^` (fix the dense `omega=...`, `R=NaN(...)`, `optProp.nin=nTIS(i)` style); split any multi-statement lines; simplify decorated `%%` banners; keep cited commented equations (e.g. the Bigio–Fantini block in `R_FD_forward.m`) fixing only their `%` spacing. Do NOT add `arguments` blocks or rename variables.

- [ ] **Step 2: Verify against baseline + checkcode**

Run:
```bash
SCRATCH="/tmp/claude-1000/-home-giles-github-gblane-dos-forward-models/c7ff4041-2dc8-466f-99e6-4254cafe2342/scratchpad"
cd /home/giles/github/gblane/dos-forward-models
matlab -batch "runVerify('$SCRATCH', 'src/analytical')"
```
Expected: all `src/analytical` regression rows `PASS` (or `SKIP(baseline errored)`), no `NEW WARNING` lines for `src/analytical` files, final line `VERIFY_OK`, exit 0. If it errors, fix the offending file and re-run before committing.

- [ ] **Step 3: Commit**

```bash
git add src/analytical
git commit -m "Conform src/analytical to MATLAB coding standard (behavior-preserving)

Whitespace, comment spacing, double-quoted strings, operator spacing, one
statement per line; public names/signatures unchanged. Verified: regression
isequaln PASS + no new checkcode ids.
Resume: next run Task 3 (src/jacobian).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: Conform `src/jacobian/` (31 files)

**Files (Modify):** all 31 `.m` files in `src/jacobian/`.

**Interfaces:** Consumes golden baseline. Produces conformed files, identical behavior.

- [ ] **Step 1: Apply the taxonomy** (same procedure as Task 2, Step 1) to every file in `src/jacobian/`. Watch the many `error('...')`, `warning('...')`, LineSpec, and property-name strings; leave `switch` `case 'DT'`/`'MC'` string labels as double-quoted (`case "DT"`) — MATLAB `switch` matches string scalars fine, but verify via regression. Keep `simTyp` default consistency: if a default is `'DT'`, make it `"DT"`.

- [ ] **Step 2: Verify against baseline + checkcode**

Run:
```bash
SCRATCH="/tmp/claude-1000/-home-giles-github-gblane-dos-forward-models/c7ff4041-2dc8-466f-99e6-4254cafe2342/scratchpad"
matlab -batch "runVerify('$SCRATCH', 'src/jacobian')"
```
Expected: all `src/jacobian` regression rows `PASS`/`SKIP`, no new warnings for `src/jacobian`, `VERIFY_OK`, exit 0.

- [ ] **Step 3: Commit** (message analogous to Task 2; resume note → "next run Task 4 (src/monte-carlo)").

---

### Task 4: Conform `src/monte-carlo/` (3 files)

**Files (Modify):** `continuousPathLen_MCadjoint.m`, `myMCXLAB_adjoint.m`, `temporalGatePathLen_MCadjoint.m`.

**Interfaces:** These are static-only (need MCXLAB/GPU + MC `adjoint` inputs). Coverage = checkcode + transformation review, reported as a coverage gap.

- [ ] **Step 1: Apply the taxonomy** to all three files.

- [ ] **Step 2: Verify (checkcode-only)**

Run:
```bash
SCRATCH="/tmp/claude-1000/-home-giles-github-gblane-dos-forward-models/c7ff4041-2dc8-466f-99e6-4254cafe2342/scratchpad"
matlab -batch "runVerify('$SCRATCH', 'src/monte-carlo')"
```
Expected: no new warnings for `src/monte-carlo`, `VERIFY_OK`, exit 0. (Regression rows for these are `SKIP` — none are in the spec registry.)

- [ ] **Step 3: Manually confirm allowed-transformation-only diff**

Run:
```bash
git diff --stat src/monte-carlo
```
Review each hunk: confirm every change is whitespace/quote/comment/operator/statement-split only — no logic change. Note in the commit that these are static-verified (no runtime coverage).

- [ ] **Step 4: Commit** (resume note → "next run Task 5 (examples/fluorescence_flow)").

---

### Task 5: Conform `examples/fluorescence_flow/` (8 files)

**Files (Modify):** `A1_MCXLAB_Homo_sen.m`, `A2_MCXLAB_Homo_sen_SD0.m`, `B1_pullPckNoi.m`, `B2_calcEta.m`, `C1_makeMap_simSig_SNRvZ.m`, `D1_realVsSimSig.m`, `D2_SNRvsNonCanNoi.m`, `run_all.m`.

**Interfaces:** Scripts (some need MCXLAB/GPU/data) → static-only: `checkcode` + parse-clean + transformation review. Plotting logic is preserved.

- [ ] **Step 1: Apply the taxonomy** to all 8 scripts. Keep all `figure/clf/plot/subplot/legend/xlabel` logic; only fix code style around it. Apply §2 path discipline lightly ONLY where already present (do not add new addpath/rmpath). `C1_makeMap_simSig_SNRvZ.m` has ~598 apostrophes — most are transposes; convert only genuine string literals, and lean on the checkcode/parse gate.

- [ ] **Step 2: Verify parse-clean + no new warnings**

Run:
```bash
SCRATCH="/tmp/claude-1000/-home-giles-github-gblane-dos-forward-models/c7ff4041-2dc8-466f-99e6-4254cafe2342/scratchpad"
matlab -batch "runVerify('$SCRATCH', 'examples/fluorescence_flow')"
```
Expected: no new warnings for these files, `VERIFY_OK`, exit 0. A parse error surfaces as a checkcode error id → gate fails.

- [ ] **Step 3: Review diff + Commit**

```bash
git diff --stat examples/fluorescence_flow
```
Confirm transformation-only. Commit (resume note → "next run Task 6 (examples/pulse_oximetry)").

---

### Task 6: Conform `examples/pulse_oximetry/` (1 file)

**Files (Modify):** `LEDspecColor_forPUB.m` (231 apostrophes, 81 trailing-WS lines, 48 no-space comments, 3 local functions with `arguments`).

**Interfaces:** Script + 3 local functions (`LEDspec_func`, `tissueOptProps_func`, `Tslab` — local copies shadowing `src/analytical`). Static-only (needs `data/MCout_avg.mat` + `makeE`); plotting preserved.

- [ ] **Step 1: Apply the taxonomy.** Convert `addpath(genpath('../../src'))`, `'linear'`, `'Color'`, LineSpecs, `'Interpreter','latex'` etc. to double quotes; fix `%nm`/`%uM`/`%1/mm` comment spacing; split `clear; home;` and `figure(N); clf;`; keep the local-function `arguments` blocks as-is. Do not touch the `.'` transposes (`SaO2_all.'`, `lam.'`).

- [ ] **Step 2: Verify parse-clean + no new warnings**

Run:
```bash
SCRATCH="/tmp/claude-1000/-home-giles-github-gblane-dos-forward-models/c7ff4041-2dc8-466f-99e6-4254cafe2342/scratchpad"
matlab -batch "runVerify('$SCRATCH', 'examples/pulse_oximetry')"
```
Expected: no new warnings, `VERIFY_OK`, exit 0.

- [ ] **Step 3: Review diff + Commit** (resume note → "next run Task 7 (my-matlab §2) + Task 8 (final sweep)").

---

### Task 7: Fix `../my-matlab/CLAUDE.md` §2 to be exemplary (separate repo)

**Files (Modify):** `../my-matlab/CLAUDE.md` (§2 code blocks only).

- [ ] **Step 1: Edit §2 code examples** so they follow §7: split `clear; home;` into two lines; change the `%% --- dependencies … ----` and `%% --- remove dependencies … ----` banners to clean `%% dependencies …` / `%% remove dependencies …`. Do not change prose meaning.

- [ ] **Step 2: Commit in the my-matlab repo**

```bash
cd /home/giles/github/gblane/my-matlab
git checkout -b claude-md-section2-conformance
git add CLAUDE.md
git commit -m "Make CLAUDE.md §2 examples conform to §7 coding standard

One statement per line (split clear; home;); plain %% section headers.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
cd /home/giles/github/gblane/dos-forward-models
```

---

### Task 8: Final full-repo verification + wrap-up

- [ ] **Step 1: Full regression + checkcode sweep**

Run:
```bash
SCRATCH="/tmp/claude-1000/-home-giles-github-gblane-dos-forward-models/c7ff4041-2dc8-466f-99e6-4254cafe2342/scratchpad"
cd /home/giles/github/gblane/dos-forward-models
matlab -batch "runVerify('$SCRATCH', 'all')"
```
Expected: every regression row `PASS`/`SKIP`, zero files with new warnings, `VERIFY_OK`, exit 0.

- [ ] **Step 2: Confirm formatting invariants repo-wide**

Run:
```bash
grep -rIlP '\t' src examples --include='*.m' || echo "no tabs"
grep -rIlP ' +$' src examples --include='*.m' || echo "no trailing whitespace"
awk 'length>120{print FILENAME": "NR}' $(find src examples -name '*.m') || echo "no long lines"
```
Expected: `no tabs`, `no trailing whitespace`, no long-line output.

- [ ] **Step 3: Update spec status + write conformance summary**

Edit `docs/superpowers/specs/2026-07-01-matlab-coding-standard-conformance-design.md` status line to `implemented`, and append a short "Verification results" note (K functions regression-PASS, static-only list, coverage gaps).

- [ ] **Step 4: Commit**

```bash
git add docs/superpowers/specs/2026-07-01-matlab-coding-standard-conformance-design.md
git commit -m "Record conformance verification results; mark spec implemented

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

- [ ] **Step 5: Report** the final state: files changed per subfolder, regression PASS count, static-only coverage gaps, and the separate `my-matlab` branch. Offer the finishing-a-development-branch skill for merge/PR.

---

## Self-Review

**Spec coverage:** §1 goal → Tasks 2–6. §2 guardrails → Global Constraints + per-task Step 1. §3 baseline gap → Task 1 Step 8 records it. §4 taxonomy → Global Constraints + Task 2 Step 1. §5 judgment calls → Global Constraints "Approved skips". §6 verification (checkcode + regression, launch check, coverage limits) → Task 1 (harness/launch), Tasks 2–6 verify steps, Task 4/5/6 static-only notes. §7 reconciliation → Global Constraints "Statements" + Task 7. §8 branch/checkpoint → per-subfolder commits with resume notes. §9 execution → workflow at handoff. §10 follow-up → Task 7. §11 out-of-scope → Global Constraints "Approved skips" + no task adds arguments/renames.

**Placeholder scan:** No TBD/TODO; every code step shows full code; every verify step shows the exact command and expected output. The `$SCRATCH` path is the concrete session scratchpad.

**Type consistency:** `regressionSpecs` returns N×4 cell `{name, fn, args, nOut}`, consumed identically by `captureBaseline` (indexes `specs{i,1..4}`) and `verifyAgainstBaseline`; `captureBaseline` returns struct `(name, ok, out, err)` consumed by `verifyAgainstBaseline`; `runCheckcode` returns `(file, ids, n)` consumed by `runVerify` via `.ids`. `runBaseline`/`runVerify` share the `outDir` baseline/checkcode `.mat` contract. Consistent.
