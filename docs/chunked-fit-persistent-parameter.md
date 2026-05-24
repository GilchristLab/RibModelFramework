# Chunked-fit pattern with a persistent parameter object

Audience: authors of adapter / runner scripts that drive AnaCoDa
(`AnaCoDa::runMCMC()`) in **chunked** mode -- many short MCMC bursts
that share state, with adaptive stop criteria (cov-diag-cv, geweke,
ar-in-band, ar-and-cov, etc.).

Applies to any model: ROC, FONSE, PA, PANSE.  Written from the PANSE
implementation, but the API points it relies on are AnaCoDa-wide.

## The problem

A naive chunked runner wraps the per-round AnaCoDa setup
(`initializeParameterObject`, `initializeModelObject`,
`initializeMCMCObject`, `runMCMC`) into a single `fitSingleRoundOf*()`
function and calls it once per chunk:

```r
repeat {
    fit <- fitSingleRoundOf*(genome, restart.input = prev.rst, samples = chunk.samples, ...)
    writeParameterObject(fit$parameter, chunk_NNN_parameter.Rdata)
    # diagnostics from fit$parameter / fit$mcmc
    # stop criterion check
    prev.rst <- last_written_rst
}
```

This re-creates the parameter and model Rcpp objects from disk on
**every** chunk.  Per-chunk fixed cost on PANSE Weinberg 0-ramp
(1781 genes, 746K positions, ~41 MB parameter object):

| Step | Cost |
|---|---|
| `initializeParameterObject(restart.file = ...)` (C++ deserialize) | 5--15s |
| `initializeModelObject(parameter, ...)` | 1--3s |
| `initializeMCMCObject(samples, thinning, ...)` | 1--2s |
| `writeParameterObject()` after run | 5--15s |
| diagnostics + small writes | 1--3s |
| **Per-chunk overhead total** | **~13--38s** |

With `chunk.samples = 100` and ~7--15 min compute per chunk, that's
5--15% of wall-clock spent on overhead.  Across hundreds of chunks
(needed for any modestly-tight stop criterion), this is the dominant
wall-clock loss above what the model itself costs.

## The pattern: hoist parameter/model out of the chunk loop

AnaCoDa's `runMCMC()` **mutates the parameter object in place** -- its
final state after `runMCMC()` returns is the natural starting state
for the next call.  This is also why `setRestartSettings()` + `.rst`
files work: the parameter's internal state at end-of-run is what
gets serialized.

The fix is to recognize that **chunks are not rounds**.  A "round"
in AnaCoDa terminology is a full conceptual segment of inference,
where you might change fix-flags, swap an NSE-rate prior, switch from
fixed-elongation to estimated, etc.  A "chunk" in an adaptive-stop
runner is just a finer-grained checkpoint **within** a single round.

So the chunked runner should look like:

```r
## ------ ONCE PER ROUND ------------------------------------------
genome    <- initializeGenomeObject(...)                  # heavy, but only once
parameter <- initializeParameterObject(genome, ...)        # heavy, but only once
parameter$fixSphi(); parameter$shareNSERate(); ...         # one-time
model     <- initializeModelObject(parameter, ...)         # one-time
model$setNSERatePriorDistribution(...)                     # one-time

## ------ PER CHUNK -----------------------------------------------
repeat {
    mcmc <- initializeMCMCObject(samples = chunk.samples, thinning,
                                 adaptive.width = adaptive.width, ...)
    mcmc$setStepsToAdapt(adaptive.steps)
    setRestartSettings(mcmc, filename = chunk.rst,
                       samples = restart.interval, write.multiple = FALSE)

    runMCMC(mcmc, genome, model, ncores)
    ## parameter is now mutated in place to end-of-chunk state.
    ## .rst was written by runMCMC for crash recovery.

    ## diagnostics
    per.aa.ar  <- capture.per.aa.ar(parameter)
    cov.diag.v <- extract.cov.diag.flat(chunk.rst)
    ## ... append to accumulators, check stop criterion ...

    ## SNAPSHOT policy: write parameter.Rdata every N chunks,
    ## not every chunk.  .rst still gives crash recovery.
    if (chunk.idx %% snapshot.every == 0) {
        writeParameterObject(parameter, chunk_NNN_parameter.Rdata)
    }
    writeMCMCObject(mcmc, chunk_NNN_mcmc.Rdata)   # small, sub-second

    if (stop_criterion_met) break
}

## ------ END OF ROUND --------------------------------------------
writeParameterObject(parameter, final_parameter.Rdata)
writeMCMCObject(mcmc, final_mcmc.Rdata)
```

Expected per-chunk overhead drops from ~13--38s to ~2--5s (just
diagnostics + the small per-chunk mcmc.Rdata write).  At
`chunk.samples = 300` and `snapshot.every = 5`, that's roughly 1--2%
of wall-clock overhead instead of 5--15%.

## Why this is safe

1. **`runMCMC` mutates parameter in place.**  This is AnaCoDa's
   standard behavior, used by every model.  After `runMCMC()`,
   `parameter` reflects end-of-chain state for the next call.
2. **Adapter state lives on `parameter`.**  CSP proposal-width
   adaptation, covariance estimates, per-codon `numAcceptForCSP`
   counters -- all stored on the parameter object's C++ side via
   `Parameter::csp_adapter` (the strategy adapter created in
   `parameter->initCSPAdapter(...)`).  Reusing parameter preserves
   tuned proposal widths across chunks, which is exactly what
   restart-from-`.rst` does -- just without the disk round-trip.
3. **Trace accumulation is bounded.**  Each `runMCMC()` call appends
   to `parameter$traces`.  At 1781 genes x 8 bytes/sample and tens of
   thousands of retained samples, this peaks in the hundreds of MB --
   well within reach for any modern fitting box.  Confirm with a
   memory check at chunk 10--20 for your specific genome size.

## Snapshot / recovery policy

The new pattern decouples three previously-coupled concerns:

| Concern | Where it lives | Frequency |
|---|---|---|
| Crash recovery (resume after kill) | `.rst` file via `setRestartSettings` | every chunk (cheap, C++ side) |
| Mid-run R-side analysis | `parameter.Rdata` via `writeParameterObject` | every N chunks (`snapshot.every`, default 5) |
| Per-chunk MCMC trace for diagnostics | `mcmc.Rdata` via `writeMCMCObject` | every chunk (small, sub-second) |

The old pattern wrote `parameter.Rdata` every chunk because each
chunk *was* a full round, so it had to.  Now that chunks are just
checkpoints within a round, the disk snapshot can be sparser.

On resume: load from the most recent `parameter.Rdata` (or the latest
`.rst` if you trust the C++ side to deserialize losslessly).
Behavior matches the old pattern's resume semantics; we just write
less often.

## What changes per model

The per-round setup function (`fitSingleRoundOf<MODEL>()`) splits
into two:

```
fitSingleRoundOf<MODEL>(genome, restart.input, samples, ...)
    ==>
initChunkSessionFor<MODEL>(genome, restart.input, ...) -> list(parameter, model)
runChunkMCMCFor<MODEL>(genome, parameter, model,
                      samples, thinning, adaptive.width, adaptive.steps,
                      ncores, divergence.iteration,
                      restart.filename, restart.interval,
                      est.csp, est.expression, est.hyper, est.mix) -> mcmc
```

The split is along the same line for every model: everything that
configures the parameter (fix flags, share.*, prior choices, init
files, init.partition.function) goes into the init function;
everything that runs MCMC goes into the run function.

### PANSE-specific

- `parameter$shareNSERate()` -- one-time
- `model$setNSERatePriorDistribution(type, lower, upper, mean,
  shape, rate)` -- one-time per round (priors don't change
  mid-round)
- `init.alpha.files`, `init.lambda.files`, `init.nserate.files`
  -- one-time (CSV-loaded init values)
- `init.partition.function` -- one-time (PANSE Z normalization)

### ROC-specific (when porting)

- `mutation.prior.*` / `selection.prior.*` arguments to
  `initializeParameterObject` -- one-time
- `parameter$fixDM()`, `parameter$fixDEta()` -- one-time per round
- `init.csp.variance` semantics (per-AA covariance) -- one-time;
  do NOT reset between chunks within a round

### FONSE-specific (when porting)

- `init.initiation.cost` -- one-time
- Stop-codon semantics (theta_i sign inversion when stop is
  included) -- handled by init.* files, one-time

## Reference implementation

PANSE: `~/Repositories/PANSE.data.analyses/s.cerevisiae/adapter.dev/lib/chunked.runner.R`
plus split helpers in `lib/local.functions.R`.  Gated behind YAML
`round.defaults.persistent.parameter` (default `false` during
validation, switch to `true` once benchmarked).

When porting to ROC / FONSE, copy the pattern:
1. Split your `fitSingleRoundOf<MODEL>()` into init + run.
2. Add the `persistent.parameter` YAML flag to your runner.
3. Move parameter/model creation before the repeat loop.
4. Inside the loop, only create the MCMC object and call runMCMC.
5. Use `snapshot.every` (default 5) for `writeParameterObject`.
6. Sanity-check trace memory at chunk 10--20 for your genome size.

## Things deliberately not addressed

- **The `writeParameterObject` cost itself.**  At ~41 MB for PANSE,
  that's 5--15s of disk I/O regardless of frequency.  Cut frequency,
  not per-call cost.  If/when needed, the lever is trace pruning
  (`parameter` carries a growing trace across chunks), not
  serialization-format changes.
- **Multi-round transitions.**  Between rounds (changing fix-flags,
  prior, etc.), the original `fitSingleRoundOf*()` path is still
  correct.  The persistent-parameter pattern applies **within** a
  single round.
- **Continuous adaptation across chunks.**  Current behavior: each
  chunk's MCMC adapts for `adaptive.steps` steps, then stops.
  Across N chunks you get N x `adaptive.steps` of adaptation.  This
  is unchanged by the refactor (adapter state is on parameter, so it
  keeps adapting), and matches the prior chunked-runner behavior
  (which also created a fresh MCMC per chunk).  If you want a one-
  shot adapt-once-then-freeze schedule, that's an orthogonal change.
