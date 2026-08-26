#!/usr/bin/env bash
#
# Launch one full-globe copdem_production.jl run, detached, for the 2026-08-26
# Conservative-vs-NearestCell comparison.
#
#     scripts/launch_copdem_bench.sh conservative
#     scripts/launch_copdem_bench.sh nearest
#
# The two runs differ ONLY in the method, the store/log names and the NUMA node
# they are pinned to. Everything else matches the phase1 reference run
# (scratch-stores/glo90-synthetic-authalic-phase1.log): -t 21 --gcthreads=4,
# workers=40, budget=2^30, schedule=:affinity, cachepolicy=:refcount,
# maskarcsec=15, resume=true, malloctrim=32MiB.
#
# NEVER pass --gcthreads=N,1: the driver's `gcguard` refuses the concurrent page
# sweeper, which segfaulted the 2026-08-21 run.
#
# Memory is capped by a transient systemd scope (a cgroup v2 memory.max), not by
# ulimit -- ulimit caps address space, not RSS, and would not stop this run.
# --heap-size-hint keeps Julia's GC working inside that cap. Note the phase1 run
# peaked at 25.2 GiB RSS, above this 24 GB cap; the run is resumable, so an
# OOM kill is recoverable by re-running this same command.
#
# RASTERDATASOURCES_PATH must point at the bench data holding
# naturalearth/ne_10m_land.shp -- the synthetic source needs it for the land
# mask, and the run dies in the first second without it.

set -euo pipefail

METHOD="${1:-}"
case "$METHOD" in
    conservative) NODE=0 ;;
    nearest)      NODE=1 ;;
    *) echo "usage: $0 conservative|nearest" >&2; exit 2 ;;
esac

RUN="glo90-synthetic-authalic-${METHOD}-20260826"
OUT=/home/asinghvi17/geo/scratch-stores
STORE="$OUT/$RUN.zarr"
LOG="$OUT/$RUN.log"
PIDFILE="$OUT/$RUN.pid"
UNIT="copdem-${METHOD}-20260826"

MEMMAX=24G
HEAPHINT=18G
THREADS=21
GCTHREADS=4

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export RASTERDATASOURCES_PATH=/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data

# `resume = true` silently skips every chunk an existing store already holds, so
# a leftover store would turn a fresh benchmark into a no-op. Refuse rather than
# measure nothing. Pass RESUME=1 to continue a run an OOM kill interrupted.
if [ -e "$STORE" ] && [ "${RESUME:-0}" != "1" ]; then
    echo "$STORE already exists; remove it (or set RESUME=1 to continue it)" >&2
    exit 1
fi
systemctl --user reset-failed "$UNIT.scope" 2>/dev/null || true

cd "$REPO"
COPDEM_METHOD="$METHOD" COPDEM_STORE="$STORE" \
setsid nohup systemd-run --user --scope \
    -p MemoryMax=$MEMMAX -p MemorySwapMax=0 --unit="$UNIT" \
    numactl --cpunodebind=$NODE --membind=$NODE \
    julia --project=benchmark -t $THREADS --gcthreads=$GCTHREADS \
        --heap-size-hint=$HEAPHINT scripts/copdem_production.jl \
    > "$LOG" 2>&1 &

# The julia PID, not systemd-run's: `numactl` execs julia in place, so the one
# process the scope's cgroup holds under this unit IS julia.
CG="/sys/fs/cgroup/user.slice/user-$(id -u).slice/user@$(id -u).service/app.slice/$UNIT.scope"
JPID=""
for _ in $(seq 1 60); do
    if [ -r "$CG/cgroup.procs" ]; then
        for p in $(cat "$CG/cgroup.procs"); do
            [ "$(cat /proc/$p/comm 2>/dev/null)" = "julia" ] && JPID=$p && break
        done
    fi
    [ -n "$JPID" ] && break
    sleep 1
done

{
    echo "pid=${JPID:-unknown}"
    echo "unit=$UNIT.scope"
    echo "cgroup=$CG"
    echo "method=$METHOD"
    echo "store=$STORE"
    echo "log=$LOG"
    echo "numanode=$NODE"
    echo "threads=$THREADS gcthreads=$GCTHREADS"
    echo "memorymax=$MEMMAX memoryswapmax=0 heap-size-hint=$HEAPHINT"
    echo "commit=$(git -C "$REPO" rev-parse HEAD)"
    echo "launched_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$PIDFILE"

echo "launched $METHOD: julia pid ${JPID:-unknown}, unit $UNIT.scope, node $NODE"
echo "  log   $LOG"
echo "  store $STORE"
echo "  pid   $PIDFILE"
