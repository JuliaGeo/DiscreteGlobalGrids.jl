# # Splitting work across tasks
#
# Every pass in this package walks `1:n` and writes each result where it
# belongs, so parallelising one is always the same move: cut `1:n` into as many
# contiguous blocks as there are tasks, and run the blocks.  Having it in one
# place is what lets a plot's `ntasks` mean the same thing everywhere — a pass
# that reached for `Threads.@threads` instead would use every thread whatever
# the plot asked for.

"""
    inparallel(f, n, ntasks; minchunk = 512)

Call `f(lo, hi)` once per contiguous block of `1:n`, on `ntasks` tasks at most.

Blocks are as even as integer division makes them, and no block is smaller than
`minchunk` — below that the spawn costs more than the work.  `ntasks = 1` runs
`f(1, n)` on the calling task, with no task spawned at all.

`f` must write only inside its own block, which is what makes the result
independent of how many tasks ran.
"""
function inparallel(f, n::Int, ntasks::Int; minchunk::Int = 512)
    n <= 0 && return nothing
    nt = clamp(ntasks, 1, max(1, cld(n, minchunk)))
    if nt == 1
        f(1, n)
        return nothing
    end
    Threads.@sync for t in 1:nt
        Threads.@spawn f(1 + div((t - 1) * n, nt), div(t * n, nt))
    end
    return nothing
end
