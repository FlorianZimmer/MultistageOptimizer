# Multistage optimizer for rockets

Calculate optimal mass for each stage given a deltaV requirement on a multistage rocket.

Finds the lowest total mass for a rocket by optimizing the weight distribution between all stages, without altering total deltaV of the rocket.

Rocket and engines are specified in a json.\
Sample files for each are in the config folder.

Mathematical background, which was also the source for this project can be found here:\
http://www.projectrho.com/public_html/rocket/multistage.php

## Config

Options can be found in config/config.json.

- `useMultiCore`: if program should utilize all available cores.
- `fullGridMode`: how the full-grid solver iterates distributions: `materialize` (default) or `streaming` (generate on the fly, avoids a giant `distributions` table; useful groundwork for GPU work but can be slower on CPU).
- `precision`: the number of discretization steps used to split total deltaV across stages (1..255).
  - Example (2 stages, precision 150): stage split candidates include `1/150` + `149/150`, `2/150` + `148/150`, ...
- `maxRAM`: maximum allowed RAM usage in bytes. Supports human-readable values like `16GB` (base 1024).
- `verbose`: debugging output (slow; only coherent in single-threaded mode).
- `enginesPath`: path to JSON where engines are defined.
- `rocketPath`: path to JSON where rocket is defined.

### Runtime budgets

- `maxCombinations`: deterministic *work* budget. If set, the program chooses the largest precision that keeps the number of evaluated distributions `<= maxCombinations` (and within `maxRAM`).
- With `zoom.enabled=true`, `maxCombinations` applies to the total work across both passes (coarse + refined). The current default split is ~70% coarse, remainder refined (the refinement radius may be reduced to stay within budget).
- `targetSeconds`: best-effort runtime target. If set (and `maxCombinations` is not set), the program runs a short calibration and derives a `maxCombinations` budget from it.
  - Note: `targetSeconds` is not deterministic across different machines (and can vary slightly even on the same machine).

### Zoom refinement (coarse-to-fine)

Enable a fast coarse search and then refine only around the best candidates at a higher precision:

- `zoom.enabled`: enable zoom mode.
- `zoom.fine_precision`: final precision (defaults to `precision`).
- `zoom.coarse_precision`: coarse pass precision (defaults to `fine_precision/2`).
- `zoom.window_coarse_steps`: neighborhood radius around the coarse winner(s), expressed in coarse steps (converted to fine units internally).
- `zoom.top_k`: refine the best K coarse candidates (reduces risk of missing the true optimum).

When zoom is enabled, the program prints `# Combinations (coarse)` and `# Combinations (refined)`. (There is no single distribution index for the refined search, so it won't print `kg at Distribution number ...`.)

### Reporting (soft warnings)

- `reporting.min_stage_dv_mps`: warn if a stage contributes less than this amount of deltaV (suggesting to remove/merge that stage). This does not affect optimization; it only prints warnings.

### Example config snippet

```json
{
  "useMultiCore": true,
  "precision": 100,
  "maxRAM": "16GB",
  "enginesPath": "config/engines.json",
  "rocketPath": "config/rocket_4stage.json",
  "maxCombinations": 5000000,
  "zoom": {
    "enabled": true,
    "fine_precision": 120,
    "coarse_precision": 60,
    "window_coarse_steps": 2,
    "top_k": 3
  },
  "reporting": {
    "min_stage_dv_mps": 200
  }
}
```

## Benchmark

Run a benchmark and emit JSON with the git head for easy comparisons:

```
MultistageOptimizer.exe --benchmark --benchmark-config config/benchmark_ci.json --benchmark-threads 1 --benchmark-iterations 3
```

You can point `--benchmark-config` at any config file (including ones that enable `zoom`, `maxCombinations`, etc.). CI uses `config/benchmark_ci.json` as a stable performance regression check.

Benchmark JSON output includes:
- `min_mass`
- `best_distribution_units` (discretized deltaV units per stage)
- `best_distribution_dv_mps` (derived deltaV per stage)

### Strategy compare wrapper

Run multiple strategies derived from the same config (full-grid baseline, zoom, budget variants) and print a combined JSON report:

```
MultistageOptimizer.exe --benchmark-compare --benchmark-config config/benchmark.json --benchmark-threads 1 --benchmark-iterations 3
```

For a human-friendly side-by-side view, use:

```
MultistageOptimizer.exe --benchmark-compare --benchmark-compare-format table --benchmark-config config/benchmark.json --benchmark-threads 1 --benchmark-iterations 3
```

**Strategies**

The compare wrapper derives strategies from the same config:

- `full_grid` (baseline): brute-force evaluate *all* distributions at the baseline precision, with `zoom`, `maxCombinations`, and `targetSeconds` removed.
  - If `zoom.enabled=true`, the baseline precision is `zoom.fine_precision` (not `precision`) so it compares against the final-resolution result; this can exceed `maxRAM` for many-stage rockets.
- `zoom`: two-pass search: full-grid at `zoom.coarse_precision`, then refine around the best `zoom.top_k` coarse candidates in a window of `zoom.window_coarse_steps` (converted to fine units) at `zoom.fine_precision`. Budgets are removed.
- `budget_full`: full-grid search with runtime budgets enabled (`maxCombinations` or `targetSeconds`), and `zoom` removed.
  - The program chooses the largest precision that fits the budget (and within `maxRAM`), which can be *higher* than `precision`.
- `budget_zoom`: zoom search with runtime budgets enabled.
  - The coarse pass uses ~70% of the budget; the refinement radius may be reduced to stay within the remaining budget.

If the baseline `full_grid` fails (usually due to RAM), comparisons are omitted and the output includes a note explaining why.

Use `--benchmark-max-seconds <sec>` to fail the run if the average total time exceeds the limit (this is what CI uses).
Benchmark runs also append a CSV row to `benchmark_results.csv` (override with `--benchmark-csv <path>`).

## Build

- Visual Studio: open `MultistageOptimizer.sln` and build.
- clang++ (quick build, no IDE):
  - `clang++ -std=c++17 -O2 -DNDEBUG -I. MultistageOptimizer.cpp -o MultistageOptimizer.exe`
  - tests: `clang++ -std=c++17 -O2 -DNDEBUG -I. -Itests tests/*.cpp -o MultistageOptimizerTests.exe`

## Releases (continuous builds)

GitHub Actions publishes rolling prereleases for every commit to `main`/`master`:

- Tag format: `b<commit-count>` (for example `b123`).
- Assets: `MultistageOptimizer-windows-x64-b<commit-count>.zip` containing `MultistageOptimizer.exe` + sample `config/`.

## Limitations

Upwards of six stages we get into realms of impossible amount of necessary RAM, with higher precisions, as the number of different distributions is calculated with n choose r where n is the precision - 1 and r is the number of stages - 1.\
To mitigate the effect, you can lower the precision.
You can try out what what number of distributions exist here: https://www.calculatorsoup.com/calculators/discretemathematics/combinations.php

## GPU acceleration (OpenCL/CUDA/ROCm)

The performance hotspot is evaluating the rocket mass for each deltaV distribution, which is embarrassingly parallel. That makes it a good *candidate* for GPU acceleration, but there’s a catch: the current full-grid solver first materializes the full `distributions` table on the CPU (`createTuple`). Offloading only the mass loop to a GPU would typically require copying `n_combinations * stages` bytes to the device (often hundreds of MB), which can wipe out most of the speedup on PCIe.

To get a meaningful GPU win, the solver generally needs to avoid building/transferring the full table and instead:
- Enumerate distributions on the device (composition generator or “index → composition” unranking) and compute mass in-kernel.
- Reduce on-device to the best mass/index (or top-k for zoom) and return only a few candidates to the CPU.

For CPU execution, materializing the distribution table can still be faster (it’s mostly a big sequential memory write, and keeps per-evaluation overhead low). The codebase includes both approaches: `optimizer::FullGridMode::Materialize` (default) and `optimizer::FullGridMode::Streaming` (on-the-fly generation, useful as a starting point for a GPU backend).

If you want to go down that route, OpenCL is the most portable single-backend option (works across NVIDIA/AMD/Intel), while CUDA/ROCm can be added as vendor-specific backends later.

Currently the program assumes that the engines always work at their vaccuum efficiency.
I tried using the sl isp of engines for the first stage only. But that always leads to a tiny first stage, because it is so "inefficient" and therefore shouldnt be very big according to the program. Because of that, the first stage engine burns for only a very short time, causing the second stage to also burn near the surface.

## KSP Mods

For using this easily and read out all necessary information in KSP, one should use:
- Kerbal Engineer Redux -> read out total mass and deltaV
- Real Fuels -> dry and wet mass of tank
- Procedural Parts -> easily and very finely adjust the size of tanks
