# `GeoMakie.GlobeAxis` draws in earth-centred coordinates, which is the one case
# where a DGGS needs no projection at all: its cells are already points on a
# sphere.  This extension recognises the axis's transform and turns it into a
# `GlobeTarget`, which reproduces what the transform would have done using
# arithmetic on the unit-sphere corner we already have.

module GeoMakieExt

import DiscreteGlobalGridsVisualization as DGGV
using DiscreteGlobalGridsVisualization: GlobeTarget, probe_ellipsoid
import GeoMakie

"""
    plot_target(tf::GeoMakie.GlobeTransform) -> GlobeTarget

The globe a `GlobeAxis` draws on.

The ellipsoid is measured out of the transform rather than assumed, so a globe
on WGS84, on a perfect sphere or in kilometres all come out at the same place
the axis's own transform would have put them.  See
[`DiscreteGlobalGridsVisualization.probe_ellipsoid`](@ref).
"""
function DGGV.plot_target(tf::GeoMakie.GlobeTransform)
    height = Float64(tf.zlevel)
    a, e2 = probe_ellipsoid(tf, height)
    return GlobeTarget(tf, a, e2, height)
end

end # module
