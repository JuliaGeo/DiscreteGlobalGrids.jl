# What the regridding methods make of the same Copernicus DEM.
#
#     julia -t auto --project=benchmark benchmark/copdem_semantics.jl \
#         res=90 box=10,11,46,47 arms=real,linear,quadratic,rim,polar \
#         methods=conservative,point
#
# `Conservative()` gives every destination cell the coverage-normalised mean of
# the source pixels it overlaps. `BarycentricPoint()` gives it a sample at its
# own centroid, interpolated between the posts around it. `NearestCell()` and
# `DirectNearest()` give it the nearest post alone. All are correct; they answer
# different questions, and this file records how far apart their answers are
# rather than asserting a bound on the distance.
#
# `scripts/copdem_production.jl` chooses between them by name, and this file
# takes the method and the missing policy from that script's own `regridmethod`
# and `regridpolicy`, so what runs here is what a production run configured the
# same way computes.
#
# ARMS, each a fresh source per method in ONE process, so the methods belong to
# one machine state. The FIRST method named is the baseline every other is
# compared against.
#
#   real       the listed GLO-90 tiles of `box` into IGeo7, on real elevations.
#              What a production box actually differs by.
#   linear     a field affine in longitude and latitude, on the scaled lattice
#              twin, with a destination FINER than the source. Q1, the trapezoid
#              and the triangle all reproduce an affine field exactly, so the
#              point arm's error against the field is roundoff; the area arm's
#              is the offset between a destination cell's centroid and the
#              centroid of the posts covering it, about half a source pixel.
#   quadratic  a field with curvature on the same lattice: no method reproduces
#              it, and the arm records what each one's error is.
#   rim        the same twin holding with a destination reaching past its edge. A
#              point stencil needing a post the holding does not have is a rim
#              and stays unmapped; an area weight over the same cell exists.
#   polar      the strip between the twin's polemost post row and the pole, which
#              no dual cell of source sample sites covers: the natural one has
#              that whole row as its corners. The point method's `poles` policy
#              decides it, and the arm runs the policy both ways against the
#              area method.
#
# WHAT IS REPORTED per arm and method: wall time of the cold sweep, the source
# chunks read, `residency`'s loads/hits/peak source bytes, peak RSS from
# `Sys.maxrss()`, the error against the analytic field where there is one, and
# against the baseline method — the fraction of cells both placed, the mean and
# maximum absolute difference there, the fraction each left missing, and how
# many cells one placed and the other did not, which for a point method against
# an area one is the rim.
#
# `benchmark/copdem_nearest.jl` owns the real source, the destination, the tile
# cache and the counting wrapper; this file includes it and reuses all of them.
#
# CONFIGURATION is `copdem_nearest.jl`'s, without `method` and `profile`, plus:
#
#   arms=a,b,c       which arms to run. Default all five.
#   methods=m,n      which methods, first one the baseline. Default
#                    conservative,point. Any name `regridmethod` takes.
#   twin=W,E,S,N     the twin's degree box. Default 10,12,49,51, which straddles
#                    the 50-degree band edge, so the transition stencils run.
#   fine=K           destination levels finer than the rule's for the twin arms.
#                    Default 1, which makes a destination cell smaller than a
#                    source post's spacing.
#
# No data is committed and no path here is absolute.

include(joinpath(@__DIR__, "copdem_nearest.jl"))

const SEMANTIC_ARMS = (:real, :linear, :quadratic, :rim, :polar)

# The scaled conformance lattice: 30 latitude intervals to the degree instead of
# 1200 or 3600, and every band ratio of a real Copernicus product. A tile of it
# is 900 posts rather than 1.4 million, so a whole arm fits in a second.
const TWIN = CD.CopernicusDEMSystem{30}()

"""
    semantics(name) -> (method, missingpolicy)

The pair `scripts/copdem_production.jl` runs under `method = name`.
"""
semantics(name) = let config = (; method = Symbol(name))
    (regridmethod(config), regridpolicy(config))
end

"""
    haspolarpolicy() -> Bool

Whether the point method in this build carries a polar policy, which it states
by taking one as a field.
"""
haspolarpolicy() = :poles in fieldnames(DGG.BarycentricPoint)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

"`copdem_nearest.jl`'s options, plus this file's own four."
function semanticoptions(args = ARGS)
    own = ("arms", "methods", "twin", "fine")
    mine = Dict{String,String}()
    rest = String[]
    for a in args
        i = findfirst('=', a)
        i === nothing && error("arguments are key=value; got $a")
        k = a[1:(i - 1)]
        k in own ? (mine[k] = a[(i + 1):end]) :
        (startswith(a, "method=") || startswith(a, "profile=") || push!(rest, a))
    end
    base = options(rest)
    arms = haskey(mine, "arms") ?
           Tuple(Symbol(s) for s in split(mine["arms"], ',')) : SEMANTIC_ARMS
    for a in arms
        a in SEMANTIC_ARMS || error("unknown arm $a; known: $(join(SEMANTIC_ARMS, ", "))")
    end
    names = haskey(mine, "methods") ? Tuple(String.(split(mine["methods"], ','))) :
            ("conservative", "point")
    length(names) >= 2 || error("methods needs a baseline and at least one other")
    foreach(semantics, names)          # reject an unknown name before any work
    twin = haskey(mine, "twin") ? Tuple(parse.(Int, split(mine["twin"], ','))) :
           (10, 12, 49, 51)
    length(twin) == 4 || error("twin is w,e,s,n in whole degrees")
    fine = haskey(mine, "fine") ? parse(Int, mine["fine"]) : 1
    return (; base..., arms, methods = names, twin, fine)
end

# ---------------------------------------------------------------------------
# Comparing two fields over the same cells
# ---------------------------------------------------------------------------

"""
    difference(base, other) -> NamedTuple

What separates two regridded fields over the same destination cells: what
fraction each placed, what fraction both placed, the spread between them where
both did, and how many cells one placed and the other did not — which for a
point method against an area one is the rim.
"""
function difference(base::AbstractVector, other::AbstractVector)
    n = length(base)
    n == length(other) || throw(DimensionMismatch(
        "the two fields cover $(n) and $(length(other)) cells"))
    both = 0
    total = 0.0
    worst = 0.0
    bmiss = 0
    omiss = 0
    bonly = 0
    oonly = 0
    @inbounds for i in 1:n
        a, p = Float64(base[i]), Float64(other[i])
        ga, gp = isfinite(a), isfinite(p)
        ga || (bmiss += 1)
        gp || (omiss += 1)
        ga && !gp && (bonly += 1)
        gp && !ga && (oonly += 1)
        (ga && gp) || continue
        both += 1
        d = abs(a - p)
        total += d
        worst = max(worst, d)
    end
    return (; n, both, bothfrac = both / n,
        baseplaced = n - bmiss, otherplaced = n - omiss,
        baseonly = bonly, otheronly = oonly,
        basemissing = bmiss / n, othermissing = omiss / n,
        meanabs = both == 0 ? NaN : total / both,
        maxabs = both == 0 ? NaN : worst)
end

"""
    fielderror(values, sites, f) -> (mean, max, n)

How far a regridded field is from the analytic field it was built from, at the
destination sample sites, over the cells it placed.
"""
function fielderror(values::AbstractVector, sites, f)
    n = 0
    total = 0.0
    worst = 0.0
    for i in eachindex(values)
        isfinite(values[i]) || continue
        n += 1
        d = abs(Float64(values[i]) - f(sites[i]))
        total += d
        worst = max(worst, d)
    end
    return (n == 0 ? NaN : total / n, n == 0 ? NaN : worst, n)
end

"Every method's field against the first one's, one block per comparison."
function reportdifferences(results, unit)
    base = results[1]
    println("\nagainst $(base.name), in $unit")
    out = Any[]
    for r in results[2:end]
        d = difference(base.out, r.out)
        @printf("  %-26s %s of %s cells (%.2f%%)\n", "$(r.name): both placed",
            commas(d.both), commas(d.n), 100 * d.bothfrac)
        @printf("  %-26s %s %.4f%%, %s %.4f%%\n", "  missing", base.name,
            100 * d.basemissing, r.name, 100 * d.othermissing)
        @printf("  %-26s %s only %s, %s only %s (%.4f%% of the %s %s placed)\n",
            "  placed by one alone", base.name, commas(d.baseonly), r.name,
            commas(d.otheronly),
            d.baseplaced == 0 ? NaN : 100 * d.baseonly / d.baseplaced,
            commas(d.baseplaced), base.name)
        @printf("  %-26s mean %.6g, max %.6g\n", "  |difference|", d.meanabs, d.maxabs)
        push!(out, (; name = r.name, d))
    end
    flush(stdout)
    return out
end

# ---------------------------------------------------------------------------
# One method over one pair
# ---------------------------------------------------------------------------

"""
    sweepmethod(name, dem, srcspace, dstspace, ndst, budget) -> NamedTuple

One method's cold sweep over every destination tile, and what it read.

The source is handed in already built, so its cache state is the caller's to
control: each method here gets one of its own, decoded from nothing.
"""
function sweepmethod(name, dem, srcspace, dstspace, ndst::Int, budget::Int)
    method, policy = semantics(name)
    println("\n$name — $(repr(method)) under $(repr(policy))")
    plan = GR.plan_regrid(dem; to = dstspace, from = srcspace, method = method,
        missingpolicy = policy, lazy = true, budget = budget)
    A = GR.regrid(dem, plan)
    # The source's own precision: rounding an affine field to Float32 costs more
    # than the interpolation that reproduces it does.
    out = Vector{eltype(A)}(undef, ndst)
    cold = phase("cold sweep", () -> sweep!(out, A, dstspace))
    st = GR.residency(A)
    reads = dem isa CountingTiles ? sum(values(dem.reads); init = 0) : nreads(dem)
    @printf("  %-26s %s chunk reads, %s loads, %s cache hits, %s peak source B\n",
        "source", commas(reads), commas(st.loads), commas(st.hits),
        commas(st.peakbytes))
    @printf("  %-26s %s of %s placed, %.2f GiB peak RSS\n", "values",
        commas(count(isfinite, out)), commas(ndst), Sys.maxrss() / 2^30)
    flush(stdout)
    return (; name, out, time = cold.time, reads, loads = st.loads, hits = st.hits,
        peakbytes = st.peakbytes)
end

# ---------------------------------------------------------------------------
# The scaled lattice twin, and analytic fields over it
# ---------------------------------------------------------------------------

"A source over a holding of twin posts that counts the reads made of it."
struct CountingPosts{C} <: DiskArrays.AbstractDiskArray{Float64,1}
    values::Vector{Float64}
    chunks::C
    reads::Threads.Atomic{Int}
end

Base.size(A::CountingPosts) = size(A.values)
DiskArrays.eachchunk(A::CountingPosts) = A.chunks
DiskArrays.haschunks(::CountingPosts) = DiskArrays.Chunked()

function DiskArrays.readblock!(A::CountingPosts, out, r::AbstractUnitRange)
    Threads.atomic_add!(A.reads, 1)
    out .= view(A.values, r)
    return out
end

nreads(A::CountingPosts) = A.reads[]

"""
    twinside(box, f) -> NamedTuple

The twin lattice's posts over `box`, held as a partial grid chunked one 1-degree
tile per chunk, carrying `f` evaluated at each post's own sample site.

The holding is a strict subset of the twin's complete level, so a stencil
reaching a post outside it is a rim — which is what the `rim` arm looks at, and
what every arm's interior is measured away from.
"""
function twinside(box, f)
    w, e, s, n = box
    tiles = sort!([Int(CD.tilecell(TWIN, lat, lon).index)
                   for lat in s:(n - 1) for lon in w:(e - 1)])
    ids = TileIds(TWIN, tiles)
    grid = DGG.PartialGrid(TWIN, 1, ids)
    space = DGG.DGGSpace(grid; chunklevel = 0)
    sites = GR.samplesites(space)
    values = [f(sites[i]) for i in 1:Int(GR.ncells(space))]
    widths = [length(GR.ownedindices(space, k)) for k in 1:Int(GR.nchunks(space))]
    chunks = DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes = widths))
    dem = CountingPosts(values, chunks, Threads.Atomic{Int}(0))
    return (; tiles, grid, space, dem, sites)
end

"""
    capdestination(sys7, south, level, chunklevel) -> NamedTuple

The IGeo7 cells of the cap north of latitude `south`, as a `PartialGrid` tiled
at `chunklevel`.

`destinationside` takes whole degrees, and the region this arm needs is a strip
one part in a few hundred of a degree deep — the gap between a lattice's
polemost post row and the pole itself.
"""
function capdestination(sys7, south::Float64, level::Int, chunklevel::Int)
    ex = Extents.Extent(X = (-180.0, 180.0), Y = (south, 90.0))
    set = DGG.query(sys7, DGG.MultiOrderCoverage(ex); level = level)
    grid = DGG.PartialGrid(DGG.CellVector(set; level = level))
    return (; grid, space = DGG.DGGSpace(grid; chunklevel = chunklevel))
end

"Longitude and latitude in degrees of a unit-sphere point."
lonlat(p) = let g = US.GeographicFromUnitSphere()(p)
    (Float64(g[1]), Float64(g[2]))
end

"An affine field in longitude and latitude, in metres."
affinefield(p) = let (lon, lat) = lonlat(p)
    1000.0 + 37.0 * lon - 61.0 * lat
end

"A field with curvature, in metres, over the same box."
curvedfield(p) = let (lon, lat) = lonlat(p)
    1000.0 + 90.0 * (lon - 11.0)^2 + 140.0 * (lat - 50.0)^2
end

"""
    polarfield(p) -> Float64

A field for the polar cap, in metres: a ramp in the plane tangent at the pole.

Longitude is not a coordinate at a pole — every meridian meets there — so a
field affine in longitude is discontinuous exactly where this arm looks, and
would report its own jump rather than the methods'. This one is affine in the
sphere's own `x` and `y`, which is affine across the cap and single-valued on
the whole sphere.
"""
polarfield(p) = 1000.0 + 7.0e5 * p[1] + 4.0e5 * p[2]

# ---------------------------------------------------------------------------
# The arms
# ---------------------------------------------------------------------------

"""
    realarm(opt, dir, sys7) -> NamedTuple

The listed GLO-90 tiles of `opt.box` into IGeo7, under every named method.
"""
function realarm(opt, dir, sys7)
    println("\n", "="^78, "\nreal — GLO-$(opt.res) $(opt.box) elevations")
    probe = sourceside(opt, dir)
    level = opt.level === nothing ? DGG.levelfor(sys7, probe.space) : opt.level
    chunklevel = opt.chunk === nothing ? max(first(DGG.levels(sys7)), level - 7) : opt.chunk
    @printf("source   %d tiles, %s posts, %.2f m spacing\n", length(probe.tiles),
        commas(DGG.ncells(probe.grid)), DGG.cellsize(probe.space))
    probe = nothing
    dst = destinationside(sys7, opt.box, level, chunklevel)
    ndst = Int(DGG.ncells(dst.grid))
    @printf("dest     IGeo7 level %d, chunk level %d, %.2f m cells, %s cells in %s tiles\n",
        level, chunklevel, DGG.cellsize(sys7, level), commas(ndst),
        commas(Int(GR.nchunks(dst.space))))
    flush(stdout)

    results = map(opt.methods) do name
        src = sourceside(opt, dir)
        sweepmethod(name, src.dem, src.space, dst.space, ndst, opt.budget)
    end
    return (; arm = :real, level, ndst, results,
        diffs = reportdifferences(results, "metres"))
end

"""
    twinarm(label, opt, sys7, f; box = opt.twin, dest = opt.twin) -> NamedTuple

An analytic field on the twin holding `box`, regridded onto the IGeo7 cells of
`dest` under every named method, and each method's own error against the field.
"""
function twinarm(label, opt, sys7, f; box = opt.twin, dest = opt.twin)
    println("\n", "="^78, "\n$label — the scaled twin over $(box), destination $(dest)")
    src = twinside(box, f)
    level = (opt.level === nothing ? DGG.levelfor(sys7, src.space) : opt.level) + opt.fine
    chunklevel = opt.chunk === nothing ? max(first(DGG.levels(sys7)), level - 3) : opt.chunk
    dst = destinationside(sys7, dest, level, chunklevel)
    ndst = Int(DGG.ncells(dst.grid))
    @printf("source   %d tiles, %s posts, %.2f m spacing in %s chunks\n",
        length(src.tiles), commas(DGG.ncells(src.grid)), DGG.cellsize(src.space),
        commas(Int(GR.nchunks(src.space))))
    @printf("dest     IGeo7 level %d, %.2f m cells, %s cells in %s tiles\n",
        level, DGG.cellsize(sys7, level), commas(ndst),
        commas(Int(GR.nchunks(dst.space))))
    flush(stdout)

    results = map(opt.methods) do name
        one = twinside(box, f)
        sweepmethod(name, one.dem, one.space, dst.space, ndst, opt.budget)
    end

    # The strict threshold against the loose one on a source with no missing
    # values: a complete point row's valid weight is its whole row, so
    # `Weighted(1)` blanks nothing `Weighted(0)` keeps.
    strict = findfirst(r -> r.name == "point", results)
    if strict !== nothing
        loose = twinside(box, f)
        looseplan = GR.plan_regrid(loose.dem; to = dst.space, from = loose.space,
            method = DGG.BarycentricPoint(), missingpolicy = DGG.Weighted(0),
            lazy = true, budget = opt.budget)
        looseA = GR.regrid(loose.dem, looseplan)
        loosevals = Vector{eltype(looseA)}(undef, ndst)
        sweep!(loosevals, looseA, dst.space)
        @printf("\n  %-26s %s under Weighted(1), %s under Weighted(0)\n",
            "point cells placed", commas(count(isfinite, results[strict].out)),
            commas(count(isfinite, loosevals)))
    end

    sites = GR.samplesites(dst.space)
    println("\nagainst the analytic field, at the destination sample sites")
    errs = map(results) do r
        m, x, n = fielderror(r.out, sites, f)
        @printf("  %-26s mean %.6g m, max %.6g m over %s placed cells\n",
            r.name, m, x, commas(n))
        (; name = r.name, mean = m, max = x, placed = n)
    end
    return (; arm = Symbol(label), level, ndst, results, errs,
        diffs = reportdifferences(results, "metres"))
end

"""
    polararm(opt, sys7) -> NamedTuple

The strip between the twin lattice's polemost post row and the pole.

No dual cell of source sample sites covers it: the natural one has that whole
row as its corners. The point method's `poles` policy decides what happens
there, and this arm runs it both ways against the area method.
"""
function polararm(opt, sys7)
    println("\n", "="^78, "\npolar — between the polemost post row and the pole")
    if !haspolarpolicy()
        println("  the point method in this build carries no polar policy, so a\n" *
                "  query poleward of the polemost post row is unmapped by\n" *
                "  construction and there is nothing to compare. Rerun once the\n" *
                "  policy is in the build.")
        return (; arm = :polar, ran = false)
    end
    # The polemost row of tiles, all the way round: the policy takes the nearest
    # site of that row, so a holding short in longitude would make its own rim
    # the thing measured instead of the policy.
    box = (-180, 180, 89, 90)
    src = twinside(box, polarfield)
    st = CD._pointstate(TWIN)
    edge = CD._sitelat(st, Int64(0))
    # The strip is a quarter of a post's latitude spacing deep, so the rule's
    # level resolves it with a cell or two. Three levels finer is 343 times the
    # cells, which is what makes a fraction a number.
    level = (opt.level === nothing ? DGG.levelfor(sys7, src.space) : opt.level) +
            opt.fine + 2
    chunklevel = opt.chunk === nothing ? max(first(DGG.levels(sys7)), level - 3) : opt.chunk
    dst = capdestination(sys7, edge, level, chunklevel)
    ndst = Int(DGG.ncells(dst.grid))
    @printf("source   %d tiles, %s posts in %s chunks; polemost row at %.6f N\n",
        length(src.tiles), commas(DGG.ncells(src.grid)),
        commas(Int(GR.nchunks(src.space))), edge)
    @printf("dest     IGeo7 level %d, %.2f m cells, %s cells north of it in %s tiles\n",
        level, DGG.cellsize(sys7, level), commas(ndst),
        commas(Int(GR.nchunks(dst.space))))
    flush(stdout)

    results = map(opt.methods) do name
        one = twinside(box, polarfield)
        sweepmethod(name, one.dem, one.space, dst.space, ndst, opt.budget)
    end

    # The same point method with the policy switched off: the difference between
    # the two point columns is the policy's whole effect.
    off = twinside(box, polarfield)
    plan = GR.plan_regrid(off.dem; to = dst.space, from = off.space,
        method = DGG.BarycentricPoint(; poles = nothing),
        missingpolicy = DGG.Weighted(1), lazy = true, budget = opt.budget)
    bareA = GR.regrid(off.dem, plan)
    bare = Vector{eltype(bareA)}(undef, ndst)
    sweep!(bare, bareA, dst.space)
    point = findfirst(r -> r.name == "point", results)
    if point !== nothing
        @printf("\n  %-26s %s of %s with poles = NearestCell(), %s with poles = nothing\n",
            "cells placed", commas(count(isfinite, results[point].out)),
            commas(ndst), commas(count(isfinite, bare)))
    end
    sites = GR.samplesites(dst.space)
    println("\nagainst the analytic field, at the destination sample sites")
    for r in results
        m, x, n = fielderror(r.out, sites, polarfield)
        @printf("  %-26s mean %.6g m, max %.6g m over %s placed cells\n",
            r.name, m, x, commas(n))
    end
    let (m, x, n) = fielderror(bare, sites, polarfield)
        @printf("  %-26s mean %.6g m, max %.6g m over %s placed cells\n",
            "point, poles = nothing", m, x, commas(n))
    end
    return (; arm = :polar, ran = true, level, ndst, results, bare,
        diffs = reportdifferences(results, "metres"))
end

# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

function run(opt = semanticoptions())
    dir = tiledir(opt)
    sys7 = DGG.IGeo7System()
    println(provenance())
    @printf("\narms=%s  methods=%s  res=GLO-%d  box=%s  twin=%s  fine=+%d  budget=%s B\n",
        join(opt.arms, ","), join(opt.methods, ","), opt.res, string(opt.box),
        string(opt.twin), opt.fine, commas(opt.budget))
    :real in opt.arms && println("tiles from $dir")

    # Compilation, on a pair of its own, so no arm below starts warm.
    warm = @timed for name in opt.methods
        w = twinside((0, 1, 0, 1), affinefield)
        lvl = DGG.levelfor(sys7, w.space)
        sp = destinationside(sys7, (0, 1, 0, 1), lvl,
            max(first(DGG.levels(sys7)), lvl - 3))
        method, policy = semantics(name)
        collect(GR.regrid(w.dem; to = sp.space, from = w.space, method = method,
            missingpolicy = policy, lazy = true, budget = opt.budget))
    end
    @printf("warm-up  %.3f s for %d method(s) (compilation; in no arm below)\n",
        warm.time, length(opt.methods))
    flush(stdout)

    out = Any[]
    for a in opt.arms
        a === :real && push!(out, realarm(opt, dir, sys7))
        a === :linear && push!(out, twinarm("linear", opt, sys7, affinefield))
        a === :quadratic && push!(out, twinarm("quadratic", opt, sys7, curvedfield))
        # The destination reaches a whole degree past the holding on every side,
        # so its edge cells need posts the holding does not have.
        a === :rim && push!(out, twinarm("rim", opt, sys7, affinefield;
            dest = (opt.twin[1] - 1, opt.twin[2] + 1, opt.twin[3] - 1, opt.twin[4] + 1)))
        a === :polar && push!(out, polararm(opt, sys7))
    end

    println("\n", "="^78, "\nsummary — each method against $(opt.methods[1])")
    @printf("%-11s %-14s %9s %10s %12s %12s %12s\n", "arm", "method", "cold s",
        "missing", "base only", "mean |diff|", "max |diff|")
    for r in out
        get(r, :ran, true) || continue
        @printf("%-11s %-14s %9.3f %9.3f%% %12s %12s %12s\n", r.arm,
            r.results[1].name, r.results[1].time,
            100 * (1 - r.diffs[1].d.baseplaced / r.diffs[1].d.n), "-", "-", "-")
        for (k, e) in enumerate(r.diffs)
            @printf("%-11s %-14s %9.3f %9.3f%% %12s %12.6g %12.6g\n", "", e.name,
                r.results[k + 1].time, 100 * e.d.othermissing,
                commas(e.d.baseonly), e.d.meanabs, e.d.maxabs)
        end
    end
    @printf("\n%d threads, peak RSS %.2f GiB\n", Threads.nthreads(),
        Sys.maxrss() / 2^30)
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
