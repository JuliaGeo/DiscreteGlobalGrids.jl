"""
Distributed CopDEM execution with a coordinator-owned graph and schedule.

Workers own their source/cache/store state and write complete, disjoint
destination chunks. Only compact completion reports return to the coordinator.
"""
module DaggerRegrid

import Dagger
import Distributed
import Printf: @sprintf
using Base.ScopedValues: @with

include(joinpath(@__DIR__, "..", "copdem_helpers.jl"))

export DaggerRegridConfig, DaggerRegridReport, dagger_smoke

const ALLOWED_FAILPOINTS = (:before_compute, :before_write, :after_write)

include("types.jl")
include("workers.jl")
include("coordinator.jl")
include("smoke.jl")

end # module DaggerRegrid
