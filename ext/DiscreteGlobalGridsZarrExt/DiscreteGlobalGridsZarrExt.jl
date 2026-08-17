# The Zarr methods behind `dggread`/`dggwrite`. This extension owns store
# access only: snapshots for the convention layer, lazy arrays for the data,
# and error enrichment with store context at the API boundary. Everything
# format-semantic (conventions, encodings, the chunked lookup) lives in the
# package's src/io layer.
module DiscreteGlobalGridsZarrExt

using DiscreteGlobalGrids
using Zarr

include("snapshot.jl")
include("read.jl")
include("write.jl")

end # module DiscreteGlobalGridsZarrExt
