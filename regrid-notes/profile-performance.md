# Profile Julia code headlessly (Tim Holy's profile-performance command)

Source: https://github.com/timholy/claude_config/blob/main/commands/profile-performance.md
Used by the P-checkpoint agents in `docs/plans/2026-08-16-regrid-tasks.md`. Analyze and report; propose optimizations, do not apply them.

## Setup

Use the MCP Julia server (`mcp__julia__julia_eval`), not `julia` via Bash — one persistent session for the whole analysis so compiled code and profile buffers survive between calls. Activate the project that owns the code as `env_path`. Profiling packages resolve from the stacked default environment, never added to the target project: `Profile` (stdlib), `FlameGraphs`, `BenchmarkTools`; optional `Cthulhu`. If missing: `using Pkg; Pkg.activate(); Pkg.add(["FlameGraphs", "BenchmarkTools"])`, then re-activate the target. Do not use ProfileView (GUI-only; FlameGraphs is its data layer).

## Step 1 — Warm up

Run the target once before measuring; the first run measures the compiler. Discard it.

## Step 2 — Baseline

```julia
using BenchmarkTools
t = @benchmark myfunc($arg1, $arg2)     # interpolate EVERY external value with $
(; t_ns = time(median(t)), bytes = memory(minimum(t)), allocs = allocs(minimum(t)),
   gc_frac = gctime(median(t)) / time(median(t)))
```

`@benchmark`, not `@btime` (prints text, returns the expression value). For long-running targets (seconds+), `@benchmark` is wrong — use `@timed` (fields `time`, `bytes`, `gctime`) for a one-shot baseline. Note whether GC is a meaningful fraction (→ Step 4) and whether the code is fast (sub-ms) or slow (→ how Step 3 samples).

## Step 3 — CPU profile

Fast code: `@bprofile myfunc($args...)` (BenchmarkTools; clears buffer itself, runs many times, GC-trial disabled). Slow code: `Profile.clear(); Profile.@profile myfunc(args...)` — `@profile` APPENDS to a global buffer, always clear before a run you analyze alone. Long jobs overflow the fixed buffer: run ≥20 s → coarsen (`Profile.init(; delay = 0.01)` or `0.05`); shorter but deep backtraces → raise `n` instead (`Profile.init(; delay=..., n = 10_000_000)`). `Profile.init()` with no args reads back current `(n, delay)`.

```julia
using Profile, BenchmarkTools, FlameGraphs
g = flamegraph()     # root Node; fully traversable, no GUI
```

Per node: `node.data.sf` (StackFrame: `.func .file .line .from_c .inlined`), `length(node.data.span)` = sample count (cost), `node.data.status` bitfield — `FlameGraphs.runtime_dispatch` (0x01, prime optimization target), `FlameGraphs.gc_event` (0x02). Iterate children with `for child in node`. Flatten and rank:

```julia
function flatten_fg(node, rows = Vector{Any}())
    total = length(node.data.span); childtotal = 0
    for c in node
        childtotal += length(c.data.span); flatten_fg(c, rows)
    end
    push!(rows, (; sf = node.data.sf, total, self = total - childtotal, status = node.data.status))
    return rows
end
rows = flatten_fg(g); sort!(rows, by = r -> -r.self)
dispatch = filter(r -> r.status & FlameGraphs.runtime_dispatch != 0, rows)
gc       = filter(r -> r.status & FlameGraphs.gc_event != 0, rows)
```

Same function appears once per call path; aggregate `self`/`total` by `(sf.func, sf.file, sf.line)` for a per-function summary. Report costliest self-time frames and, separately, costliest runtime-dispatch frames. On Julia 1.12+ use `@profile_walltime` for I/O-bound or `@spawn`-heavy code (CPU profiler under-counts waiting).

## Step 4 — Allocation profile

```julia
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=0.1 myfunc(arg1, arg2)   # 1.0 exact/slow, 0.01 long runs
allocs = Profile.Allocs.fetch().allocs   # Vector of Alloc: .type .size .stacktrace (leaf first)
```

Aggregate bytes by `.type` and by `.stacktrace[1]`, sort descending. Julia 1.10 reports some as `UnknownType` (expected); filter `CorruptType`/`BufferType` sentinels if noisy.

## Step 5 — Type instabilities

For the costliest dispatch frames: `@inferred f(args...)` (binary check), then `@code_warntype f(args...)` — look for `::Any`, `::Union{...}`, `Core.Box` (closure capturing a reassigned variable). Follow chains by running `@code_warntype` on the unstable callee. Cthulhu only if recursive callsite analysis is genuinely needed (version-fragile).

## Step 6 — Report

Baseline numbers; top CPU self-cost frames with `file:line`; top runtime-dispatch frames + instabilities behind them; top allocation types and sites; concrete prioritized proposals, each tied to evidence, biggest-win first. After an accepted change, re-benchmark in the same session and `judge(median(after), median(before))` to confirm it beats noise.

Agent's edge over the human workflow: don't eyeball the flame graph — exhaustively traverse and aggregate (cost per function, ranked dispatch sites, bytes per type). Lean into that.
