# Architecture — Non-interactive analysis interface (PC / inverse-normal combination)

## Context
The exported function `adaptGMCP_PC()` performs multi-look graphical multiple testing with closed testing and inverse-normal p-value combination across looks. Today it is implemented as an interactive console workflow (via `readline()`), which blocks:

- Batch processing of multiple analysis cases
- Writing fully scripted, end-to-end analysis pipelines
- Unit testing against deterministic inputs

## Goals
- Preserve backward compatibility: **do not modify** existing exported functions.
- Provide a fully non-interactive analysis interface.
- Implement a stateless pipeline: each step takes an explicit state object + new inputs and returns updated state.

## Non-goals
- Refactoring the legacy interactive `adaptGMCP_PC()` implementation.
- Changing the core mathematics.
- Implementing the CER method in this change (can be done later using the same pattern).

## Design overview
### API surface
Three new exported functions implement the workflow:

- `SetupAnalysis_PC()`
  - Creates a "PCAnalysisState" object.
  - Computes design boundaries (rpact), intersection weights (`genWeights()`), and inverse-normal weights.

- `AnalyzeLook_PC(state, p_raw, ...)`
  - Advances the analysis by exactly one look.
  - Optional (look > 1) pre-analysis updates:
    - `selection` (drop/retain hypotheses)
    - `new_weights` + `new_G` (strategy update)
    - `new_correlation` (correlation update)

- `PlotAnalysisGraph(state, stage)`
  - Plots the graph from the current state or from historical snapshots.

Additionally, S3 methods exist for convenience:
- `print.PCAnalysisState()` — prints design + per-look results so far
- `plot.PCAnalysisState()` — plots the current graph

### Stateless workflow
All mutable information is carried in the returned object. There is no hidden global state.

Recommended usage:
1. `state <- SetupAnalysis_PC(...)`
2. For each look `k`:
   - `state <- AnalyzeLook_PC(state, p_raw = ..., selection = ..., new_weights = ..., new_G = ..., new_correlation = ...)`

## State object
Class: `PCAnalysisState` (S3)

Key fields:
- `mcpObj`: internal list mirroring the legacy `mcpObj` structure used by `PerLookMCPAnalysis()`.
- `setup_mcpObj`: snapshot of the setup state to support plotting the initial graph.
- `thresholds`: stage-wise p-value thresholds.
- `mvtnorm_algo`: algorithm object from `chooseMVTAlgo()`.
- `completed_looks`: integer.
- `trial_completed`: logical — TRUE when no active hypotheses remain.
- `look_history`: list of per-look snapshots (stores the `mcpObj` after each analysis and the inputs used).

## Termination behavior
The pipeline stops when further analysis is impossible.

Hard stops:
- `trial_completed == TRUE` (all hypotheses rejected or dropped)
- `completed_looks == LastLook` (all planned looks completed)

Unlike the interactive workflow (which has user-driven continuation), the non-interactive workflow does not prompt; users control the flow by deciding whether to call `AnalyzeLook_PC()` again.

## Reuse of existing core computations
The non-interactive interface deliberately reuses existing internal functions without modification:
- `PerLookMCPAnalysis()`
- `addNAPvalue()`
- `modifyIntersectWeights()`
- `genWeights()`
- `chooseMVTAlgo()`
- `plotGraph()`

This limits risk and ensures numerical parity with `adaptGMCP_PC()`.

## Testing approach
The legacy function is interactive and cannot be called directly from automated tests.

Strategy:
- Add **test scaffolds** for the new API using deterministic inputs.
- Gate tests behind an option (skipped by default), and leave TODO placeholders for expected outputs.
- A developer can run `adaptGMCP_PC()` manually for each test scenario and paste expected results into the tests, then enable execution.

## Files
- R/pc_analysis_state.R
- R/pc_analysis_helpers.R
- R/pc_analysis_api.R
- tests/testthat/test-pc_analysis_api.R
- internalData/PC_Analysis_NonInteractive_Example.R
