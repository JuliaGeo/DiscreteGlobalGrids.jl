# ---------------------------------------------------------------------------
# The budget traversal
#
# Same target preparation, prunes and exact predicates as `_coverage_visit!`;
# what changes is the SCHEDULE and what stops it.
#
# The covering is {contained cells} ∪ {crossing members}. Contained cells are
# never refined. Crossing members refine COARSEST FIRST, ties broken by curve
# order within the level, so the schedule is a function of the inputs alone.
# The queue is two vectors and a pass counter: each pass sorts both streams on
# `(level, key)` and advances one level; members refine before phantoms. A
# cell kept whole is never opened — a member beneath it would be emitted under
# its own ancestor.
#
# PHANTOMS carry the descent through cells that miss the target but pass both
# `node_extent` prunes: under non-congruent refinement a cell can meet the
# target while its parent misses it. Phantoms are free; the meeting cells
# found beneath them are ENTRANTS and cost one whole slot each. On a congruent
# system every meeting cell has a meeting parent, so the member stream reaches
# everything and phantoms stay empty.
#
# PENDING cells are members none of whose children meet: the target lies in
# their overhang annulus and the phantom stream reaches the cells covering it.
# A pended cell leaves the covering but keeps its slot as a CLAIM until
#   * REDEEMED — at the margin, by an entrant that provably stands in for it
#     (see `_claim_match`);
#   * COMPLETED — its neighbourhood searched to the end, every meeting cell
#     admitted: the slot is released to the stream;
#   * RESOLVED — by `_budget_resolve!`: dropped where the search covered its
#     share, handed back where the evidence stands on its cell.
# Only non-congruent systems pend.
#
# A replacement that does not fit leaves its cell kept whole. An entrant that
# can pay neither a fresh slot nor a claim is REFUSED and recorded; after a
# refusal, every branch no unsettled claim keeps alive dies as ABANDONED —
# territory never revisited, which keeps later refusals local and lets settled
# claims release soundly. Refused and abandoned cells together bound
# everything the walk left uncovered: the EVIDENCE the end phase settles
# claims against.
# ---------------------------------------------------------------------------

function _multi_order_budget(sys::AbstractHierarchicalGridSystem, target_value,
        maxcells::Int, maxlevel::Int)
    maxcells >= 1 || throw(ArgumentError(
        "maxcells must be at least 1, got $maxcells"))
    maxlevel in levels(sys) || throw(ArgumentError(
        "maxlevel $maxlevel is outside $(typeof(sys))'s levels $(levels(sys))"))
    target = _query_target(target_value)
    ID = cellindextype(sys)
    top = first(levels(sys))
    grids = [levelgrid(sys, l) for l in top:maxlevel]

    # Congruent systems need no phantoms: every meeting cell has a meeting parent.
    congruent = has_congruent_refinement(sys)

    contained = ID[]        # proven `Within`: never refined, never requeued
    stalled = ID[]          # crossing, and kept whole: the budget said no
    current = ID[]          # crossing members, at the level being refined
    phantoms = ID[]         # non-members that passed the prunes, same level
    pending = ID[]          # members no child of which meets: settled at the end
    settled = Bool[]        # per pending cell: no longer holds a budget slot
    covered = Bool[]        # per pending cell: share proven covered
    evidenced = Bool[]      # per pending cell: a refusal recorded on its cell
    refused = ID[]          # meeting cells the budget could not pay for
    abandoned = ID[]        # branches a refusal killed: unsearched territory
    entrant_refused = false
    for c in rootcells(sys)
        # The prune is shared: `_budget_admit!` assumes it ran, and rejection
        # feeds the phantom stream.
        _budget_reaches(sys, target, c) || continue
        _budget_admit!(contained, current, sys, target, c, grids, top, maxlevel) && continue
        congruent || push!(phantoms, c)
    end
    # The budget's running claim: contained + crossing + stalled + unsettled
    # pending. A pended cell holds its claim until it is settled — redeemed,
    # released on completion, or resolved at the end — so the end phase always
    # has room for what is still uncovered. Phantoms are outside it.
    total = length(contained) + length(current)
    schedule_key(c) = (level(c), _budget_key(grids[level(c)-top+1], c))
    cellcap(c) = cell_cap(grids[level(c)-top+1], c)
    # What a redeeming entrant must contain whole (see `_claim_match`).
    pieces = congruent ? nothing : _target_pieces(target_value, target)

    kids_in = ID[]
    kids_out = ID[]
    kids_miss = ID[]
    for _ in top:(maxlevel-1)
        (isempty(current) && isempty(phantoms)) && break
        sort!(current; by=schedule_key)
        sort!(phantoms; by=schedule_key)
        next = ID[]
        next_phantoms = ID[]
        for c in current
            # Unreachable while both streams enter at the roots and advance one
            # level per pass; keeps `grids` in range if entry ever deepens.
            if level(c) >= maxlevel
                push!(stalled, c)
                continue
            end
            empty!(kids_in)
            empty!(kids_out)
            empty!(kids_miss)
            for child in children(sys, c)
                _budget_reaches(sys, target, child) || continue
                _budget_admit!(kids_in, kids_out, sys, target, child, grids,
                    top, maxlevel) && continue
                congruent || push!(kids_miss, child)
            end
            k = length(kids_in) + length(kids_out)
            if k == 0 && !congruent
                # No child meets: the target lies in this cell's overhang
                # annulus, and the phantom stream reaches the cells that cover
                # it. The claim is held only where handing the cell back could
                # ever be the right answer — commensurate with a piece it
                # guards; a cell that dwarfs every piece it touches frees its
                # slot to the stream, and the refused entrants stand in for
                # its share at the end.
                push!(pending, c)
                held = _claim_worth_holding(cellcap(c),
                    ncells(grids[level(c)-top+1]), pieces)
                push!(settled, !held)
                push!(covered, false)
                push!(evidenced, false)
                held || (total -= 1)
                append!(next_phantoms, kids_miss)
            elseif k == 0 || total + k - 1 > maxcells
                # Kept whole; the subtree stays closed so no member lands under
                # a member. Congruent `k == 0` lands here too: children that
                # tile their parent cannot lose the target between them — the
                # cell meets it on a boundary sliver — and dropping it would
                # break the exactness this mode promises where refinement is.
                push!(stalled, c)
            else
                append!(contained, kids_in)
                append!(next, kids_out)
                total += k - 1
                append!(next_phantoms, kids_miss)
            end
        end
        # The unsettled claims, deepest first: what marginal entrants redeem,
        # and — after a refusal — the only territory still worth searching.
        live = [i for i in eachindex(pending) if !covered[i]]
        sort!(live; by=i -> (-level(pending[i]), schedule_key(pending[i])[2]))
        livecaps = [cellcap(pending[i]) for i in live]
        for p in phantoms
            # After a refusal, a branch matters only where a claim's share is
            # still unproven. Around an EVIDENCED claim it matters only for
            # what lies beyond that claim's own cell: one refusal there
            # already justifies the hand-back, which covers the cell whole,
            # and searching inside it refines nothing the end phase can use.
            if entrant_refused
                ext = node_extent(sys, p)
                keep = false
                for k in eachindex(live)
                    i = live[k]
                    covered[i] && continue
                    Extents.intersects(livecaps[k], ext) || continue
                    if !evidenced[i] || !_cell_contains_cap(
                            cell_cap(grids[level(p)-top+1], p),
                            grids[level(pending[i])-top+1], pending[i])
                        keep = true
                        break
                    end
                end
                if !keep
                    push!(abandoned, p)
                    continue
                end
            end
            # Same `grids`-range guard as the member loop.
            level(p) >= maxlevel && continue
            for child in children(sys, p)
                _budget_reaches(sys, target, child) || continue
                empty!(kids_in)
                empty!(kids_out)
                if !_budget_admit!(kids_in, kids_out, sys, target, child, grids,
                        top, maxlevel)
                    push!(next_phantoms, child)
                    continue
                end
                # A meeting cell whose parent misses replaces nothing, so it
                # costs a whole cell: a fresh slot, or — at the margin — the
                # slot of a claim it provably stands in for.
                if total + 1 <= maxcells
                    total += 1
                elseif (m = _claim_match(grids[level(child)-top+1], child,
                        pieces, live, livecaps, settled, total, maxcells)) != 0
                    settled[m] = true
                    covered[m] = true
                else
                    entrant_refused = true
                    push!(refused, child)
                    rcap = cellcap(child)
                    for k in eachindex(live)
                        i = live[k]
                        (covered[i] || evidenced[i]) && continue
                        _cell_meets_cap(rcap, grids[level(pending[i])-top+1],
                            pending[i], false) && (evidenced[i] = true)
                    end
                    continue
                end
                append!(contained, kids_in)
                append!(next, kids_out)
            end
        end
        # A claim's share is proven covered by COMPLETION when nothing can
        # change its covering any more: no live phantom still reaches it, and
        # no refusal was recorded against it — nor can one be later, since
        # any branch able to reach it would still be live now. Everything
        # found there was admitted, so the share is covered; a slot still
        # held is surplus, released for the stream to spend.
        if !congruent && !isempty(pending)
            for i in eachindex(pending)
                covered[i] && continue
                cap = cellcap(pending[i])
                any(p -> Extents.intersects(cap, node_extent(sys, p)),
                    next_phantoms) && continue
                any(a -> Extents.intersects(cap, node_extent(sys, a)),
                    Iterators.flatten((abandoned, refused))) && continue
                covered[i] = true
                if !settled[i]
                    settled[i] = true
                    total -= 1
                end
            end
        end
        current = next
        phantoms = next_phantoms
    end
    append!(stalled, current)       # members only: what `maxlevel` left unrefined
    _budget_resolve!(stalled, contained, sys, pending, settled, covered,
        grids, top, maxcells, abandoned, refused, pieces)

    cells = vcat(contained, stalled)
    flags = falses(length(cells))
    flags[1:length(contained)] .= true
    # The reference level is the deepest level reached; an empty set reports the top.
    reference = isempty(cells) ? top : maximum(level, cells)
    return _sorted_cell_set(sys, cells, flags, reference)
end

# `node_extent` bounds the whole subtree, so a cell failing these prunes
# shelters nothing that meets the target. The same test `_coverage_visit!`
# descends on.
function _budget_reaches(sys, target, c)
    extent = node_extent(sys, c)
    Extents.intersects(target.cap, extent) || return false
    return !_subtree_outside(target, extent)
end

# Classify one candidate: `Intersects`, then `Within` — the predicates
# `_coverage_visit!` uses, in the same order; at `maxlevel`, `Within` is skipped
# (see `iscontained`). Callers run `_budget_reaches` first: a rejected cell is a
# phantom candidate, and both decisions share that one prune.
function _budget_admit!(contained, crossing, sys, target, c, grids, top::Int,
        maxlevel::Int)
    lc = level(c)
    grid = grids[lc-top+1]
    _matches(DE9IM.Intersects(nothing), target, grid, c) || return false
    if lc < maxlevel && _matches(DE9IM.Within(nothing), target, grid, c)
        push!(contained, c)
    else
        push!(crossing, c)
    end
    return true
end

# One (bounding cap, target) pair per connected piece of the target — a
# MultiPolygon's parts; any other target is one piece, itself. Parts prepare
# on first use: most walks never test containment.
mutable struct _PieceSlot{G}
    const geom::G
    prepared::Any
    points::Any
end

function _target_pieces(target_value, target)
    if GI.isgeometry(target_value) &&
       GI.trait(target_value) isa GI.MultiPolygonTrait
        return [(points_cap([query_point(p) for p in GI.getpoint(g)]),
                 _PieceSlot(g, nothing, nothing)) for g in GI.getpolygon(target_value)]
    end
    target isa GeometryTarget || return [(target.cap, target)]
    return [(target.cap, _PieceSlot(target_value, target, nothing))]
end

# A held claim may come back whole at the end — only ever right for a cell
# commensurate with a piece it touches; a giant handed back over a speck is
# the annulus pathology, and its slot serves the stream better. Yardstick:
# the piece cap's spherical area against the cell's, tuned on the annulus and
# refusal-isolation laws (a cap already overshoots its piece, so 1 carries
# the slack a covering cell needs).
const _CLAIM_PROPORTION = 1.0

_cap_fraction(cap) = (1 - cos(min(Float64(cap.radius), Float64(pi)))) / 2

function _claim_worth_holding(cellcap, ncells_level, pieces)
    yardstick = 0.0
    for (cap, _) in pieces
        Extents.intersects(cellcap, cap) || continue
        yardstick = max(yardstick, _cap_fraction(cap))
    end
    return 1 / ncells_level <= _CLAIM_PROPORTION * yardstick
end

# Whether one whole piece of the target lies inside the cell — the only
# containment the engine answers exactly; a piece kind without an exact test
# certifies nothing. The slot rejects on the first vertex outside the cell's
# cap and prepares its geometry only when a candidate survives that.
_piece_in_cell(piece::GeometryTarget, grid, c, cellcap) =
    _matches(DE9IM.Contains(nothing), piece, grid, c)
function _piece_in_cell(slot::_PieceSlot, grid, c, cellcap)
    if slot.points === nothing
        g = GI.isgeometry(slot.geom) ? slot.geom : _extent_target(slot.geom)
        slot.points = [query_point(p) for p in GI.getpoint(g)]
    end
    for p in slot.points
        cap_contains(cellcap, p) || return false
    end
    slot.prepared === nothing && (slot.prepared = _query_target(slot.geom))
    return _piece_in_cell(slot.prepared::GeometryTarget, grid, c, cellcap)
end
_piece_in_cell(piece, grid, c, cellcap) = false

# The deepest unsettled claim the entrant provably stands in for, or 0: the
# entrant's CELL must contain every piece the claim's cap touches. A cap test
# is not enough — a cap swallows more than its cell, and a claim redeemed on
# a sliver loses the rest of its share. `live` comes deepest first, so the
# first hit is the smallest cell. Redemption is a transfer, not a purchase:
# refused while the covering itself is over budget (a seed above `maxcells`).
function _claim_match(grid, child, pieces, live, livecaps, settled, total::Int,
        maxcells::Int)
    total <= maxcells || return 0
    isempty(live) && return 0
    childcap = cell_cap(grid, child)
    holds = zeros(UInt8, length(pieces))    # 0 unknown, 1 contains, 2 does not
    for k in eachindex(live)
        settled[live[k]] && continue
        ok = false
        for (j, (cap, piece)) in enumerate(pieces)
            Extents.intersects(livecaps[k], cap) || continue
            if holds[j] == 0x00
                holds[j] = _piece_in_cell(piece, grid, child, childcap) ? 0x01 : 0x02
            end
            ok = holds[j] == 0x01
            ok || break
        end
        ok && return live[k]
    end
    return 0
end

# Hands back what the walk still owes, finest cells first, inside `room`.
#
# Everything uncovered lies inside a refused entrant or under a culled
# branch: those are the evidence. Candidates are the refused cells (each its
# own evidence) and the pending cells not proven covered, which need evidence
# ON THEIR CELL — not their cap, which reaches territory the footprint never
# held: a giant handed back on far-away evidence is the annulus pathology.
#
# One pass, deepest first: emitting a refused cell settles its evidence, so a
# coarser pending cell comes back only for evidence nothing finer explained,
# and of two nested candidates the smaller blocks the coarser.
#
# STARVED pieces come before everything: a piece no member stands on has only
# this pass to be covered at all, and room spent first on refused slivers of
# a well-covered piece loses the island whole.
#
# A covering is never empty: with nothing else kept, the shallowest unproven
# cell — the one whose footprint holds the most target — is the answer.
function _budget_resolve!(stalled, contained, sys, pending, settled, covered,
        grids, top::Int, maxcells::Int, abandoned, refused, pieces)
    unresolved = [pending[i] for i in eachindex(pending) if !covered[i]]
    (isempty(unresolved) && isempty(refused)) && return nothing
    room = maxcells - (length(contained) + length(stalled))
    ev = vcat(refused, abandoned)
    evcaps = [cell_cap(grids[level(a)-top+1], a) for a in ev]
    members = Iterators.flatten((contained, stalled))
    evlive = [j > length(refused) ||
              !_refused_covered(sys, ev[j], evcaps[j], members, grids, top)
              for j in eachindex(ev)]
    if any(evlive) && room > 0
        # Held claims only: a freed cell was freed because handing it back
        # could never be proportionate (see `_claim_worth_holding`).
        cands = [(pending[i], 0) for i in eachindex(pending) if !settled[i]]
        append!(cands, [(refused[j], j) for j in eachindex(refused) if evlive[j]])
        # A NEGLECTED piece — far fewer members than the dominant one — is a
        # component the budget passed over, and losing it whole is the
        # multipart pathology; its candidates come first, coarsest first,
        # covering the most with the least. Everything else is tightening,
        # finest first.
        mcount = zeros(Int, length(pieces))
        for m in members
            for j in eachindex(pieces)
                _cell_meets_cap(pieces[j][1], grids[level(m)-top+1], m, false) &&
                    (mcount[j] += 1)
            end
        end
        mmax = maximum(mcount; init=0)
        rescues(c) = any(j -> 2 * mcount[j] < mmax && Extents.intersects(
                cell_cap(grids[level(c)-top+1], c), pieces[j][1]),
            eachindex(pieces))
        sort!(cands; by=((c, _),) -> begin
            r = rescues(c)
            (!r, r ? Int(level(c)) : -Int(level(c)),
             _budget_key(grids[level(c)-top+1], c))
        end)
        for (c, j) in cands
            room > 0 || break
            (_ancestor_of_any(sys, c, contained) ||
             _ancestor_of_any(sys, c, stalled) ||
             _descends_from_any(sys, c, contained) ||
             _descends_from_any(sys, c, stalled)) && continue
            if j == 0
                grid = grids[level(c)-top+1]
                any(k -> evlive[k] && _cell_meets_cap(evcaps[k], grid, c, false),
                    eachindex(ev)) || continue
            end
            push!(stalled, c)
            j == 0 || (evlive[j] = false)
            room -= 1
        end
    end
    (isempty(contained) && isempty(stalled) && !isempty(unresolved)) &&
        push!(stalled, argmin(level, unresolved))
    return nothing
end

# Emitting both a cell and its ancestor claims the same leaves twice; consumers
# read a set as disjoint subtrees.
function _ancestor_of_any(sys, c, members)
    lc = level(c)
    return any(m -> level(m) > lc && ancestor(sys, m, lc) == c, members)
end

# Whether the cap lies entirely inside the cell: centre in the cell, and every
# boundary arc at least the radius away. The converse of `_cell_meets_cap`,
# with the same exactness for convex caps.
function _cell_contains_cap(cap, grid, c)
    ring, n = open_ring(cell_boundary(grid, c))
    point_in_cell(ring, cap.point) === true || return false
    threshold = cos(min(Float64(pi), Float64(cap.radius)))
    for i in 1:n
        a = ring[i]
        b = ring[i == n ? 1 : i+1]
        arc = BoundaryArc(USPoint(a[1], a[2], a[3]), USPoint(b[1], b[2], b[3]))
        _arc_cos_distance(arc, cap.point) > threshold && return false
    end
    return true
end

# A refused cell some single member's footprint holds whole is already
# covered: not a candidate, and not evidence either.
function _refused_covered(sys, e, ecap, members, grids, top)
    for m in members
        mcap = cell_cap(grids[level(m)-top+1], m)
        US.spherical_distance(mcap.point, ecap.point) <= mcap.radius || continue
        _cell_contains_cap(ecap, grids[level(m)-top+1], m) && return true
    end
    return false
end

function _descends_from_any(sys, c, members)
    lc = level(c)
    return any(m -> level(m) < lc && ancestor(sys, c, level(m)) == m, members)
end

# The tie-break inside one level. `globalindex` is the level's own order, which
# is curve order on every system here and is a bijection, so no two cells of a
# level ever tie and the schedule has nothing left to decide. A missing index
# would silently alias two cells onto one key and break determinism — so it is
# an error, not a zero.
function _budget_key(grid, c)::Int
    pos = globalindex(grid, c)
    pos === nothing && throw(ArgumentError(
        "$(typeof(grid)) has no index for the cell $c it just produced. The " *
        "budget schedule orders each level by `globalindex` and needs it total: " *
        "without it two cells share a key and the traversal stops being " *
        "deterministic"))
    return pos
end
