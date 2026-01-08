# Open issues (remaining feedback)

This file tracks the **remaining items** from the review feedback ("DELIVERABLE_FOR_CODEX_CLI" plan).

Each issue includes:
- **Problem**: what's missing / why it matters
- **Fix**: concrete implementation checklist
- **Done when**: acceptance criteria + suggested validation commands

Priority legend:
- **P0** must-have (correctness / safety / baseline tooling)
- **P1** should-have (performance + KSP alignment)
- **P2** nice-to-have / future feature

## OI-001 (P0) Sanitizers + warnings build mode (Plan Step 1)
- [ ] **Status**: open

**Problem**
- There is no repeatable "warnings + ASan/UBSan" build/run path, so UB/range regressions can slip in unnoticed.

**Fix**
- Add **one** repeatable way to build and run with warnings + sanitizers:
  - Option A (recommended): add `scripts/build_asan.ps1` that builds:
    - `optimizer_asan.exe` from `MultistageOptimizer.cpp`
    - `tests/test_distributions_asan.exe`, `tests/test_calcMass_asan.exe`, `tests/test_json_import_asan.exe`
  - Option B: add a `Sanitizers` (or similar) configuration in `MultistageOptimizer.vcxproj` using clang-cl/LLVM.
- Use (or equivalent):
  - `-std=c++17 -O1 -g -fno-omit-frame-pointer`
  - `-fsanitize=address,undefined`
  - `-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion`
- Ensure non-interactive execution (use `--no-prompt` or `--benchmark`).

**Touches**
- new: `scripts/build_asan.ps1` (or `MultistageOptimizer.vcxproj`)

**Done when**
- `optimizer_asan.exe --benchmark --benchmark-config config/benchmark.json --benchmark-iterations 1 --benchmark-warmup 0` runs with **no ASan/UBSan reports**.
- All test executables run with no ASan/UBSan reports.

## OI-002 (P1) Performance: remove `massDistributions` array (Plan Step 13)
- [ ] **Status**: open

**Problem**
- The solver currently:
  - allocates `massDistributions[nCombinations]`
  - fills it
  - scans it again to find the minimum
- This costs memory and adds an extra full pass over `nCombinations`.

**Fix**
- Switch to "per-thread best" + reduction:
  - Define a result type, e.g.:
    - `struct Best { double mass; std::uint64_t index; std::vector<unsigned char> distro; };`
  - Add a function that computes best-in-range without storing all masses, e.g.:
    - `Best distBest(bool verbose, uint64_t begin, uint64_t end, const Rocket&, const unsigned char* distributions, int precision, ...)`
  - Spawn threads so each thread writes exactly one `Best` (unique slot, no mutex).
  - Reduce to global best with deterministic tie-break:
    - lower mass wins
    - if equal (within epsilon), lower index wins
  - Remove:
    - the `massDistributions` allocation
    - the post-pass min scan
  - Update RAM checks to reflect actual allocations:
    - keep `distributionsBytes`
    - remove `resultsBytes` for `massDistributions`
    - add small overhead for `best_per_thread` + the final `bestDistro`

**Touches**
- `math/massCalc.hpp`
- `MultistageOptimizer.cpp`

**Done when**
- Output best distribution + mass matches a reference run for the same inputs.
- Benchmark shows reduced memory usage and improved runtime (especially `min_seconds`, often `total_seconds` too).

## OI-003 (P1) KSP alignment: `g0`, `ispMode`, optional loss model (Plan Step 14)
- [ ] **Status**: open

**Problem**
- No config-driven KSP-like Isp handling:
  - always uses `engine.isp_vac` (no `sl_first_stage`)
  - no consistent `g0` control (KSP uses `g0 = 9.80665` for converting Isp seconds -> exhaust velocity)
  - no dV loss model

**Fix**
- Extend config schema (`config/config.json`) with:
  - `"g0": 9.80665`
  - `"ispMode": "vac" | "sl_first_stage"` (optional later: `"blend"`)
  - optional `"deltaVLoss": 0` (simple additive loss model)
  - optional `"validateDistributions": true|false` (so distribution validation can run without `verbose`)
- Update `MultistageOptimizer.cpp` to read these fields and pass them into mass evaluation.
- Update `math/massCalc.hpp` to select stage Isp based on `ispMode`:
  - `vac`: always `engine.isp_vac`
  - `sl_first_stage`: use `engine.isp_sl` for the **bottom** stage only, vacuum for others
- Implement `deltaVLoss` consistently (choose one and document it in code/README):
  - Option A: treat as extra required dV -> `deltaVmission += deltaVLoss`
  - Option B: treat as unavailable dV -> `deltaVmission -= deltaVLoss` (must not go negative)
- Add/extend a minimal test/invariant to ensure toggling `ispMode` changes results in the expected direction.

**Touches**
- `config/config.json`
- `MultistageOptimizer.cpp`
- `math/massCalc.hpp`
- tests: extend `tests/test_calcMass.cpp` (or add a new test)

**Done when**
- Running with `ispMode=sl_first_stage` produces a higher (worse) optimal mass than `ispMode=vac` for typical rockets.
- Results remain deterministic across repeated runs for each mode.

## OI-004 (P2) Asparagus + side staging foundation (burn segment model) (Plan Step 15)
- [ ] **Status**: open

**Problem**
- Current solver is strictly serial stages. Side staging / asparagus requires:
  - parallel burns (multiple stacks burning simultaneously)
  - jettison events mid-ascent
  - effective Isp for multiple engines, which requires thrust data

**Fix**
- Extend engine data and types (required for meaningful multi-engine Isp):
  - Add `thrust_sl` / `thrust_vac` to `config/engines.json`
  - Add those fields to `rocket_definition/Engine.hpp`
  - Import them in `read_json/read_json.hpp`
- Add a separate (isolated) parallel-staging model:
  - `rocket_definition/Stack.hpp` (core/booster definition)
  - `rocket_definition/BurnSegment.hpp` (which stacks burn, crossfeed rules, what gets jettisoned)
  - `optimizer/ParallelStagingSolver.hpp` (or similar)
- Implement segment effective Isp (thrust-weighted harmonic mean):
  - `Isp_eff = F_total / sum_i(F_i / Isp_i)`
- Implement in this order:
  1) side staging without crossfeed (boosters + core burn together, drop boosters)
  2) asparagus crossfeed (outer feeds inner/core), modeled as repeated segments
  3) only after model is stable: integrate optimization (segment dV allocation + prop sizing)

**Touches**
- `config/engines.json`
- `rocket_definition/Engine.hpp`
- `read_json/read_json.hpp`
- new: `rocket_definition/Stack.hpp`
- new: `rocket_definition/BurnSegment.hpp`
- new: `optimizer/ParallelStagingSolver.hpp`

**Done when**
- A standalone solver can compute mass/dV for a fixed parallel-staging design and produces reasonable results.
- Existing serial-stage optimization remains unchanged and continues to work.

## OI-005 (P2) Repo hygiene: ignore generated binaries
- [ ] **Status**: open

**Problem**
- Local builds create untracked binaries (e.g. `optimizer.exe`, `tests/*.exe`).

**Fix**
- Update `.gitignore` to ignore:
  - `optimizer.exe`
  - `tests/*.exe`
  - (optional) `optimizer_asan.exe`, `tests/*_asan.exe` once OI-001 exists

**Touches**
- `.gitignore`

**Done when**
- `git status` is clean after building + running tests (excluding intentional source edits).

## OI-006 (P2) Style/portability: remove `using namespace std;` from `MultistageOptimizer.cpp`
- [ ] **Status**: open

**Problem**
- `using namespace std;` in a large TU increases name-collision risk and makes refactors noisier.

**Fix**
- Remove `using namespace std;` from `MultistageOptimizer.cpp`.
- Update code to use `std::` explicitly (or use narrow `using std::string;` etc. in local scopes).

**Touches**
- `MultistageOptimizer.cpp`

**Done when**
- Builds successfully; benchmark outputs are unchanged.

