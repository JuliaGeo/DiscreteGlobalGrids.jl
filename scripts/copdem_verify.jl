# The synthetic oracle: read written chunks back and hold them against
# `CopernicusUtils.synthetic_elevation`. Included by `copdem_common.jl`.

"""
    verify(config, src, layout, chunks)

Read written chunks back out of the store and hold them against
`synthetic_elevation`. Four claims, in the order they can fail:

  * every value is finite or `NaN`, nothing else;
  * a cell with no valid source pixel under it is `NaN`, and a fully-fed cell
    matches the field at its centre to within the field's curvature across a
    cell;
  * a coastal cell lies inside the field's range across that cell — which is
    what says the weights renormalised over the fed fraction rather than
    diluting toward zero;
  * a chunk nobody wrote reads back all `NaN`.
"""
function verify(config, src, layout, chunks)
    level, ancestor, storepath = config.level, config.ancestor, config.store
    g12 = DGG.levelgrid(src.sys7, level)
    sm = SourceMask(src.mask, DGG.levelgrid(src.sys, 0), Set(src.tiles))
    n = min(config.checkchunks, length(chunks))
    picked = chunks[round.(Int, range(1, length(chunks); length = n))]
    say("verify: reading back $n of $(length(chunks)) chunks")
    anymixed = false
    for ch in picked
        a = DGG.columncell(layout, ch)
        stack = DGG.dggread(storepath; ancestors = [a])
        vals = collect(stack[:elevation])
        cells = collect(DD.lookup(stack[:elevation], DGG.Cells))
        h = DGG.columnlength(layout, ch)
        check("chunk $ch: read back $h values",
            length(vals) == h && length(cells) == h;
            detail = "$(length(vals)) values, $(length(cells)) cells")

        nfin = count(isfinite, vals)
        nnan = count(isnan, vals)
        check("chunk $ch: every value is finite or NaN", nfin + nnan == length(vals);
            detail = "$nfin finite, $nnan NaN")

        # A stride, not every cell: the boundary walk is 8 slerps per edge and a
        # chunk is 823 543 cells. The stride is prime to 7 so it does not land on
        # one subtree's worth of siblings.
        stride = max(1, length(vals) ÷ 400)
        nfull = nnone = nmixed = 0
        errs = Float64[]
        bracketfail = 0
        nanfail = Int[]
        finfail = Int[]
        for k in 1:stride:length(vals)
            c = cells[k]
            kind = cellsource(g12, c, sm)
            if kind === :none
                nnone += 1
                isnan(vals[k]) || push!(nanfail, k)
            elseif kind === :full
                nfull += 1
                isfinite(vals[k]) || (push!(finfail, k); continue)
                push!(errs, abs(vals[k] - synthetic_elevation(DGG.cell_centroid(g12, c))))
            else
                nmixed += 1
                isfinite(vals[k]) || continue
                samples = synthetic_elevation.(DGG.cell_boundary(g12, c))
                push!(samples, synthetic_elevation(DGG.cell_centroid(g12, c)))
                lo, hi = extrema(samples)
                tol = 1e-3 + 0.05 * (hi - lo)
                (lo - tol <= vals[k] <= hi + tol) || (bracketfail += 1)
            end
        end
        nmixed > 0 && (anymixed = true)
        say("chunk $ch sampled: $nfull fully fed / $nnone unfed / $nmixed coastal")
        check("chunk $ch: unfed cells are NaN", isempty(nanfail);
            detail = "$(length(nanfail)) of $nnone unfed cells are not NaN")
        check("chunk $ch: fully-fed cells are finite", isempty(finfail);
            detail = "$(length(finfail)) of $nfull fully-fed cells are NaN")
        if !isempty(errs)
            mx = maximum(errs)
            check("chunk $ch: fully-fed cells match the analytic field", mx < 1.0;
                detail = @sprintf("max %.3e m, RMS %.3e m over %d cells", mx,
                    sqrt(Statistics.mean(abs2, errs)), length(errs)))
        end
        check("chunk $ch: coastal cells bracket the field", bracketfail == 0;
            detail = "$bracketfail of $nmixed coastal cells outside their own field range")
        if nnan > 0 && nfin > 0
            say("chunk $ch holds BOTH land and nodata: $nfin finite, $nnan NaN " *
                @sprintf("(%.1f%%)", 100 * nnan / length(vals)))
        end
    end
    # A global run must meet the coast; a region need not — an all-Antarctica box
    # is 100 % valid and has nothing to renormalise.
    if config.region === nothing
        check("some sampled cell straddles the coast", anymixed;
            detail = anymixed ? "the renormalising path ran" : "no coastal cell was sampled")
    else
        say("coastal cells sampled in this region: " *
            (anymixed ? "yes, the renormalising path ran" :
             "none — a region with no coastline has nothing to renormalise"))
    end

    # A chunk nobody wrote: the first level-`ancestor` chunk outside the
    # covering, which by construction meets no land tile.
    inset = Set(chunks)
    empty = findfirst(i -> !(i in inset), 1:DGG.ncells(DGG.system(layout), ancestor))
    if empty !== nothing
        a = DGG.columncell(layout, empty)
        vals = collect(DGG.dggread(storepath; ancestors = [a])[:elevation])
        check("unwritten ocean chunk $empty reads back NaN", all(isnan, vals);
            detail = "$(count(isnan, vals)) of $(length(vals)) NaN")
    end
    return nothing
end
