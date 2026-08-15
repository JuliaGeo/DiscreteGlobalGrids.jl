# Aggregation over the hierarchy: `aggregate` reduces leaf data to one fixed
# coarser level; `coarsen` merges complete sibling groups within a tolerance,
# bottom-up, into a `MultiOrderVector` + values — the adaptive-mesh
# constructor. `expand` is the inverse presentation (its DimArray methods live
# in `src/dimensionaldata.jl`; the function is born here).
#
# Contract: docs/design/moc-storage.md §2.

function aggregate end
function coarsen end
function expand end
