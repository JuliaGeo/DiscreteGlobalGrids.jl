# `MultiOrderVector` — the mixed-level cell container: cells at different
# refinement levels, pairwise-disjoint subtrees, sorted and indexed by their
# descendant-range intervals at a reference level. The storage form of a MOC,
# where `MultiOrderCellSet` is the query form.
#
# Contract: docs/design/moc-storage.md §1.

# The complement of a multi-order container over its system's whole sphere.
# Born here as a zero-method function so the export binds; methods below.
function complement end
