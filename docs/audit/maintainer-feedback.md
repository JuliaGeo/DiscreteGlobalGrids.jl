# Maintainer review decisions

These notes record feedback after the audit. They do not authorize package implementation changes.
The original sweep reports remain historical evidence; these decisions supersede their recommendations where they conflict.

- Remove `systems()` and the registry in a later implementation pass. Check internal callers when planning removal.
- Keep S2 and Copernicus DEM out of the high-level analysis tutorial. Availability alone does not justify equal prominence.
- Do not teach direct `PartialGrid` construction as an ordinary user workflow. Drop API-002, the proposed convenience constructor.
- The remainder of `sweep-core.md` is generally accepted, subject to these decisions.
- Keep API-003, but focus it on a public, consistent `cell_polygon` operation returning GeoInterface polygons. `cell_boundary` remains the vertex-level operation. Expanding `cell_polygons` is not the agreed objective.
- Explain API gaps with a concrete task, current behavior, desired behavior, and a short code example. Label proposed syntax as unimplemented.

See [the revised API gaps](api-gaps.md) for the active proposals. Implementation remains deferred while review continues.
