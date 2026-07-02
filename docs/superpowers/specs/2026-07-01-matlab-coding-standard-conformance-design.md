# MATLAB coding-standard conformance pass — design

Date: 2026-07-01
Branch: `matlab-coding-standard`
Status: **implemented** (2026-07-02) — all 57 files conformed and verified

## 0. Verification results (2026-07-02)

All 57 `.m` files (`src/` 48 + `examples/` 9) conformed and committed. Final
full-repo gate: **VERIFY_OK** — 24/24 regression-guarded functions bit-identical
(`isequaln`) to the pristine baseline, **0 regressions, 0 new checkcode ids**.
Formatting invariants repo-wide: 0 tabs, 0 trailing-whitespace lines, 0 lines
> 120 chars. Static-only files (2L/Scat/Grad, `makeSenMaps*`, `makeS`,
`monte-carlo/*`, all example scripts) additionally verified behavior-preserving
by a normalized code-skeleton comparison (`verify/skeleton.py`: comments/
whitespace stripped, quote delimiters unified → identical HEAD vs conformed).

Notable decisions found during implementation:
- `makeE` chromophore codes (`'OD'`/`'O'`/`'D'`/`'ODWL'`) and MCX `cfg` flags
  (`'flux'`/`'p'`/`'cone'`/`'pencil'`) stay **char** (indexed/parsed
  char-by-char); bracket char-concatenations `['a' x 'b']` / `warning([...])`
  stay char (double quotes → string array = behavior change). The regression
  harness caught an over-conversion of `makeE('ODWL')`.
- Example scripts conformed via `verify/conform_quotes.py` (deterministic,
  transpose/bracket/continuation-aware). Section-banner decoration in example
  scripts was left as-authored (behavior-neutral; deviates from §5 "simplify
  banners", kept for consistency across `examples/` and minimal churn).
- Pre-existing checkcode warnings on static-only files (NANSUM/NANMEAN, AGROW,
  CAXIS, NASGU, UNRCH) were left as-is — the gate is "no *new* warnings," and
  these fixes (deprecated-fn renames, preallocation, dead-code removal) are
  behavior-adjacent and not runtime-verifiable here. Recommended follow-ups.
- Task 7 (`../my-matlab/CLAUDE.md` §2) was moot: that standard was restructured
  into `../my-matlab/AGENTS.md`, whose examples already conform (no
  `clear; home;`, no decorated banners).
- `verify/` (MATLAB harness + `conform_quotes.py` + `skeleton.py`) is branch
  tooling — keep or drop before merge.

## 1. Goal

Bring every MATLAB source file in this repository into conformance with the
MATLAB coding standard defined in `../my-matlab/CLAUDE.md` §7 (which defers to the
MathWorks MATLAB Coding Guidelines). The pass covers all 57 `.m` files under
`src/` (analytical, jacobian, monte-carlo) and `examples/` (fluorescence_flow,
pulse_oximetry).

The enforceable bar is: **Code-Analyzer-clean, no new `checkcode` warnings (fewer
where safely achievable), plus the formatting / string / comment best-practices** —
achieved **without changing behavior, public function names, signatures, or
plotting output**.

## 2. Guardrails (hard constraints)

- **No public renames.** Function file names and signatures are preserved exactly.
  Many (`R_FD_forward`, `PHI_TD_forward`, `Tslab`, `n2A`, `LEDspec_func`,
  `continuousReflectance`, `makeSenMaps`, …) are the canonical API that the sibling
  repos (`dos-inverse-models`, `dual-slope-toolkit`) and the examples call by name.
  Renaming would break downstream callers (this repo is their canonical home — see
  `../my-matlab/CLAUDE.md` §5).
- **No behavior change.** Every edit is behavior-preserving and proven per file by a
  `checkcode` pass plus an output-equality regression check — asserted with
  evidence, not by inspection.
- **No plotting rework.** Example scripts keep their existing `figure/clf/plot/
  subplot` logic; they are not converted to the `saveFigure` house style (§2 of the
  standard). Only the code style *around* plotting is fixed.
- **Hazard avoidance.** Transpose operators (`'`, `.'`) and in-string tokens (`%d`,
  `%f` inside `sprintf`/`fprintf`) are never altered by the quote or comment passes.
  A blind `sed`/regex sweep would corrupt these — edits are semantic-aware and
  regression-verified.

## 3. Baseline gap (measured 2026-07-01)

- Tabs: **0** (indentation is already spaces).
- Lines > 120 chars: **0** (line length already conforms).
- Trailing-whitespace lines: **~510**.
- Comments missing the space after `%` (mostly unit comments `%mm`, `%1/mm`,
  `%rad/sec`): **~608**.
- Apostrophes: **~1599** — but this count mixes real single-quoted string literals
  with transpose operators, so the true string-conversion count is lower and must
  be found semantically, not by count.
- `arguments` blocks: present in most `src/jacobian`, `src/monte-carlo`, and newer
  `src/analytical` functions; absent in older analytical functions that use
  `nargin`-based defaults.

## 4. Change taxonomy (applied per file, safest first)

1. **Whitespace** — strip trailing whitespace; confirm 4-space indent and no tabs;
   remove leading/trailing blank lines in each file; single blank lines between
   logical blocks and before local functions.
2. **Comments** — exactly one space after `%` (comments only; never touch `%%`
   section markers, `%{`/`%}` block comments, or `%`-tokens inside string literals).
3. **Strings** — convert single-quoted **string literals** to `"double-quoted"`,
   skipping transposes and any call site that genuinely requires a char vector.
   Correctness is confirmed by regression, not by the regex alone.
4. **Operator spacing** — spaces around `=`, relational, and logical operators; no
   spaces around `:`/`*`/`/`/`^`, inside brackets, in `Name=Value`, or after unary
   `+`/`-`/`~` (fixes the older `omega=...`, `R=NaN(...)` dense style).
5. **Statements** — **strict one statement per line**: split `clear; home;`,
   `figure(101); clf;`, `subplot(2,1,1); hold on;`, etc. (This tightens the earlier
   "keep idiomatic pairs" idea; see §7 reconciliation below.)
6. **Code Analyzer** — resolve remaining safe `checkcode` warnings (unused
   variables, `~` for dropped outputs, etc.) without altering behavior.

## 5. Judgment calls (approved defaults)

- **Dead / reference code** (e.g. the commented Bigio–Fantini equation in
  `R_FD_forward`): **keep** — these are cited alternate derivations, not dead code.
  Only their comment spacing is normalized.
- **`%% ###…###` / `%% ---…---` decorated section banners**: **simplify** to clean
  `%% Title`.
- **`nargin` → `arguments` block conversion**: **skip** — it changes the input
  contract (shape/type coercion, validation) and risks behavior change. Existing
  `arguments` blocks are kept; no new ones are added.
- **Variable renaming** (`Mfrac_all`, `SaO2_all`, `m_max` → camelCase): **skip** —
  not Analyzer-detectable, high churn with no test net, and domain acronyms (`SaO2`)
  read worse in camelCase. Local variables are only touched when a `checkcode`
  warning requires it.

## 6. Verification strategy (MATLAB is installed at `/usr/local/bin/matlab`)

1. **Confirm MATLAB launches** (license/headless) before any bulk edit; if it does
   not, downgrade to `checkcode`-only verification and record the reduced coverage.
2. **Regression harness (behavior-preservation).** For every function callable with
   computable default inputs, capture its output on fixed inputs **before** the
   refactor, then **after**, and assert `isequal` (or `isapprox` with a tight
   tolerance for floating-point) — this is the primary behavior-preserving proof.
3. **`checkcode` gate.** Run `checkcode` on every file before and after; require no
   new warnings; record warnings resolved.
4. **Coverage limits (stated, not hidden).** Example scripts that depend on MCXLAB /
   GPU (`examples/fluorescence_flow/A1`, `A2`, and the pulse-oximetry Monte-Carlo
   path) cannot be fully executed here. They receive `checkcode` + parse-clean
   verification only, and this reduced coverage is reported per file.

## 7. Reconciliation with the standard document

The standard's own `../my-matlab/CLAUDE.md` §2 currently shows `clear; home;` on one
line and `%% --- … ---` decorated banners, which conflict with §7's
"one statement per line" and clean-section-header rules. We resolve the tension in
favor of §7 (strict one-statement-per-line, clean banners) in this repo, and record
a **follow-up side task** (see §10) to make the standard document's §2 examples
exemplary so the two no longer disagree.

## 8. Branch, checkpointing, and resumability (standard §6)

- Work on branch `matlab-coding-standard`.
- **Commit the verification harness first**, before the first bulk edit, so an
  interruption during editing loses nothing.
- Commit **per subfolder** — `src/analytical`, `src/jacobian`, `src/monte-carlo`,
  `examples/fluorescence_flow`, `examples/pulse_oximetry` — each with a short resume
  note (what is done, what is next, the exact next command).
- Each commit is independently verified (its files pass `checkcode` + regression)
  so a resumed session restarts at the first unconformed subfolder.

## 9. Execution model

Once the implementation plan is approved, run the pass as a verification-gated
workflow: a per-file pipeline of `snapshot outputs → conform edits → checkcode +
regression verify`, each file independent, with failures reported rather than
hidden. Subfolder commits follow each verified batch.

## 10. Follow-up side task (separate repo)

Update `../my-matlab/CLAUDE.md` §2 so its example code is itself exemplary of §7:
split `clear; home;` into one statement per line and simplify the `%% --- … ---`
banners to clean `%% Title`. This is a distinct edit/commit in the `my-matlab`
repository, tracked separately from this repo's conformance branch.

## 11. Out of scope

- Renaming public functions or changing signatures.
- Converting example plotting to the `saveFigure` house style.
- Adding a `CLAUDE.md` governance file to this repo.
- Adding `arguments` blocks to functions that currently use `nargin` defaults.
- Sweeping variable renames.
