isdefined(@__MODULE__, :DaggerTalkFigures) ||
    include(joinpath(@__DIR__, "11-dagger-theme.jl"))

using .DGGSTalkFigures, .DaggerTalkFigures
using Bonito, DiscreteGlobalGridsVisualization, Makie, WGLMakie
import SparseArrays

const TIMESERIES_OUTPUT = html_path("13-regrid-timeseries.html")

# The export contains the small synthetic cube in JavaScript. Once the page is
# written, the browser interpolates the adjacent slices itself; no Julia process
# is involved in either playback or scrubbing.
const TIMESERIES_STYLE = Styles(
    CSS("html, body", "width" => "100%", "height" => "100%", "margin" => "0",
        "display" => "grid", "place-items" => "center", "overflow" => "hidden",
        "background" => "#f8f9fa"),
    CSS(".regrid-timeseries", "position" => "relative", "width" => "960px",
        "height" => "540px", "overflow" => "hidden", "background" => "#ffffff"),
    # Reserve a real bottom band for the control.  Letting the slider float over
    # the canvas made its white background cut through the lower edge of both
    # globes at slide size.
    # Three independently sized canvases are deliberately composed by CSS rather
    # than by one responsive Makie GridLayout.  A globe is a square camera
    # viewport, and letting the source, plan, and destination negotiate widths
    # made their centres drift when Slidev resized the embedded page.
    CSS(".regrid-timeseries__panel", "position" => "absolute", "overflow" => "hidden"),
    CSS(".regrid-timeseries__source", "left" => "34px", "top" => "42px",
        "width" => "300px", "height" => "328px"),
    CSS(".regrid-timeseries__plan", "left" => "405px", "top" => "116px",
        "width" => "150px", "height" => "220px"),
    CSS(".regrid-timeseries__destination", "left" => "626px", "top" => "42px",
        "width" => "300px", "height" => "328px"),
    CSS(".regrid-timeseries__heading", "position" => "absolute", "left" => "0",
        "top" => "0", "width" => "100%", "height" => "20px",
        "display" => "flex", "align-items" => "center", "justify-content" => "center",
        "font-family" => "Inter, sans-serif", "font-size" => "10px",
        "font-weight" => "600", "letter-spacing" => "0.03em", "white-space" => "nowrap",
        "color" => "#212529"),
    CSS(".regrid-timeseries__plan .regrid-timeseries__heading", "color" => "#2c7a1e"),
    CSS(".regrid-timeseries__canvas", "position" => "absolute", "left" => "0",
        "top" => "22px", "width" => "100%", "height" => "calc(100% - 22px)"),
    CSS(".regrid-timeseries__timeline", "position" => "absolute", "z-index" => "4",
        "left" => "122px", "right" => "36px", "bottom" => "17px", "height" => "24px",
        "display" => "flex", "align-items" => "center", "gap" => "10px",
        "padding" => "0 10px", "border" => "1px solid #dee2e6",
        "border-radius" => "3px", "background" => "rgba(255,255,255,0.90)",
        "box-shadow" => "0 1px 4px rgba(0,0,0,0.09)"),
    CSS(".regrid-timeseries__time", "font-family" => "Inter, sans-serif",
        "font-size" => "10px", "font-weight" => "600", "letter-spacing" => "0.08em",
        "white-space" => "nowrap", "color" => "#6c757d"),
    CSS(".regrid-timeseries__slider", "width" => "100%", "accent-color" => "#389826",
        "opacity" => "0.82"),
)

"""
    regrid_timeseries_app(; nt = 12)

Build a static WGLMakie/Bonito page in which a 24 × 12 raster globe and its
IGeo7 regridded globe advance together. The sparse plan is immutable and only
the two color buffers are updated in browser-side JavaScript.
"""
function regrid_timeseries_app(; nt = 12)
    fixture = regrid_fixture(; nt)
    weights = fixture.plan.block.weights
    rows, cols, vals = SparseArrays.findnz(weights)
    # A travelling analytic field makes the time dimension legible: one pass
    # moves its crests through 360° of longitude. Applying the already-built
    # sparse weights here keeps the destination cube an actual regridding of
    # the source rather than a second, independently drawn field.
    source_cube = hcat([
        vec([1000sind(3(lon - 360 * (k - 1) / nt)) * cosd(2lat) + 100
            for lon in fixture.lon, lat in fixture.lat])
        for k in 1:nt
    ]...)
    # `block.weights` carries area-scaled contributions. The plan normalizes
    # each destination row when applying it, so preserve that same weighted
    # average here while keeping the exported operator visibly unchanged.
    row_weight = vec(sum(weights; dims = 2))
    destination_cube = (weights * source_cube) ./ row_weight
    source_cells = longlat_cells(length(fixture.lat))

    # `vec` makes colors a browser-updatable buffer. The operator is rendered
    # once, while the two geographic color observables change entirely in JS.
    source_colors = Observable(copy(source_cube[:, 1]))
    destination_colors = Observable(copy(destination_cube[:, 1]))
    source_frames = [collect(source_cube[:, k]) for k in 1:nt]
    destination_frames = [collect(destination_cube[:, k]) for k in 1:nt]

    with_theme(JG_THEME) do
        # Fixed figure dimensions are mirrored exactly around the compact
        # operator.  The HTML below gives these canvases their fixed locations
        # inside the 960 × 540 export.
        source_fig = Figure(; size = (300, 306), figure_padding = 0,
            backgroundcolor = JG.paper)
        source = globe_axis(source_fig[1, 1]; camera_longlat = (20, 18),
            camera_altitude = 2.50)
        plot_cells!(source, source_cells; color = source_colors,
            colormap = FIELD_COLORMAP, colorrange = FIELD_RANGE,
            strokecolor = (JG.ink, 0.18), strokewidth = 0.24)
        coastlines!(source)

        plan_fig = Figure(; size = (150, 198), figure_padding = (3, 2, 6, 3),
            backgroundcolor = JG.paper)
        plan = Axis(plan_fig[1, 1]; xlabel = "source", ylabel = "IGeo7", aspect = 0.62)
        scatter!(plan, cols, rows; color = vals, colormap = FIELD_COLORMAP,
            colorrange = extrema(vals), markersize = 1.55)
        xlims!(plan, 1, size(weights, 2)); ylims!(plan, size(weights, 1), 1)
        plan.xticks = ([1, size(weights, 2)], ["1", string(size(weights, 2))])
        plan.yticks = ([1, size(weights, 1)], ["1", string(size(weights, 1))])

        destination_fig = Figure(; size = (300, 306), figure_padding = 0,
            backgroundcolor = JG.paper)
        destination = globe_axis(destination_fig[1, 1]; camera_longlat = (20, 18),
            camera_altitude = 2.50)
        dggpoly!(destination, fixture.grid; color = destination_colors,
            colormap = FIELD_COLORMAP, colorrange = FIELD_RANGE,
            strokecolor = (JG.ink, 0.34), strokewidth = 0.52)
        coastlines!(destination)

        # Plot insertion otherwise lets the two unrelated meshes pick separate
        # final zooms.  Commit the same camera after all layers are present.
        for axis in (source, destination)
            axis.center[] = false
            cameracontrols(axis.scene).settings.center[] = false
            Makie.update_cam!(axis; longlat = (20, 18), altitude = 2.50,
                fov = 45.0)
        end

        slider = DOM.input(; class = "regrid-timeseries__slider", type = "range",
            min = "0", max = "1000", step = "1", value = "0",
            ariaLabel = "Time interpolation")
        time_readout = DOM.span("SLICE 01 / $(lpad(nt, 2, '0'))";
            class = "regrid-timeseries__time")
        root = DOM.div(
            DOM.div(
                DOM.div("SOURCE RASTER · 24 × 12 × TIME";
                    class = "regrid-timeseries__heading"),
                DOM.div(WGLMakie.WithConfig(source_fig; resize_to = :parent);
                    class = "regrid-timeseries__canvas");
                class = "regrid-timeseries__panel regrid-timeseries__source"),
            DOM.div(
                DOM.div("SPARSE PLAN · ONCE"; class = "regrid-timeseries__heading"),
                DOM.div(WGLMakie.WithConfig(plan_fig; resize_to = :parent);
                    class = "regrid-timeseries__canvas");
                class = "regrid-timeseries__panel regrid-timeseries__plan"),
            DOM.div(
                DOM.div("IGEO7 DESTINATION · 492 CELLS × TIME";
                    class = "regrid-timeseries__heading"),
                DOM.div(WGLMakie.WithConfig(destination_fig; resize_to = :parent);
                    class = "regrid-timeseries__canvas");
                class = "regrid-timeseries__panel regrid-timeseries__destination"),
            DOM.div(time_readout, slider; class = "regrid-timeseries__timeline");
            class = "regrid-timeseries",
        )

        App() do session
            Bonito.onload(session, root, js"""async (root) => {
                const slider = root.querySelector('.regrid-timeseries__slider');
                const readout = root.querySelector('.regrid-timeseries__time');
                const sourceFrames = $(source_frames);
                const destinationFrames = $(destination_frames);
                const sliceCount = sourceFrames.length;

                if (!slider || !readout || !sliceCount) {
                    console.error('regrid timeseries: offline controls unavailable');
                    return;
                }

                let phase = 0.0;
                let previous = performance.now();
                let scrubbing = false;
                const secondsPerRotation = 9.0;

                const blend = (a, b, t) => a.map((v, i) => (1.0 - t) * v + t * b[i]);

                // A Julia Observable reaches WGLMakie through Julia's compute
                // graph.  In an exported page there is no Julia process to run
                // that graph, so updating the source Observable alone only
                // changes its JS value.  Locate the two numeric vertex-color
                // buffers in WGLMakie's public cache and update those GPU-bound
                // attributes directly instead.
                function closest_indices(attribute, values) {
                    const indices = new Uint16Array(attribute.array.length);
                    let maxError = 0.0;
                    for (let i = 0; i < attribute.array.length; i++) {
                        const value = attribute.array[i];
                        let best = 0, error = Infinity;
                        for (let j = 0; j < values.length; j++) {
                            const candidate = Math.abs(value - values[j]);
                            if (candidate < error) {
                                best = j;
                                error = candidate;
                            }
                        }
                        indices[i] = best;
                        maxError = Math.max(maxError, error);
                    }
                    return { indices, maxError };
                }

                function field_buffer(frames) {
                    let winner = null;
                    for (const mesh of Object.values(window.WGL?.plot_cache || {})) {
                        const attribute = mesh?.geometry?.attributes?.vertex_color;
                        if (!attribute || attribute.itemSize !== 1 ||
                            attribute.array.length < frames[0].length) continue;
                        const match = closest_indices(attribute, frames[0]);
                        if (!winner || match.maxError < winner.maxError) {
                            winner = { attribute, ...match };
                        }
                    }
                    return winner?.maxError < 0.02 ? winner : null;
                }

                async function wait_for_field_buffers() {
                    for (let attempt = 0; attempt < 120; attempt++) {
                        const source = field_buffer(sourceFrames);
                        const destination = field_buffer(destinationFrames);
                        if (source && destination && source.attribute !== destination.attribute) {
                            return { source, destination };
                        }
                        await new Promise(requestAnimationFrame);
                    }
                    throw new Error('regrid timeseries: WGL color buffers unavailable');
                }

                const fields = await wait_for_field_buffers();
                // Buffer discovery may take several animation frames on first
                // load; start playback from the moment the canvas is ready.
                previous = performance.now();
                const screens = Array.from(root.querySelectorAll('canvas'))
                    .map(canvas => canvas.wglmakie_screen)
                    .filter(Boolean);
                function update_field(field, values) {
                    const { array } = field.attribute;
                    for (let i = 0; i < array.length; i++) array[i] = values[field.indices[i]];
                    field.attribute.needsUpdate = true;
                }

                function setTime(u) {
                    const wrapped = Math.max(0.0, Math.min(0.999999, u));
                    const progress = wrapped * sliceCount;
                    const first = Math.floor(progress);
                    const alpha = progress - first;
                    const second = (first + 1) % sliceCount;
                    const nextSource = blend(sourceFrames[first], sourceFrames[second], alpha);
                    const nextDestination = blend(destinationFrames[first], destinationFrames[second], alpha);
                    update_field(fields.source, nextSource);
                    update_field(fields.destination, nextDestination);
                    // `needsUpdate` uploads the changed attributes on the next
                    // Three.js render; static WGLMakie otherwise has no Julia
                    // scene event to request that frame.
                    for (const screen of screens) {
                        if (screen?.root_scene) window.WGL.render_scene(screen.root_scene);
                    }
                    slider.value = String(Math.round(wrapped * 1000));
                    readout.textContent = `SLICE ${String(first + 1).padStart(2, '0')} / ${String(sliceCount).padStart(2, '0')}`;
                    // Keep the active value observable to visual-regression
                    // checks without adding another display element to slide.
                    root.dataset.frame = `${first}:${alpha.toFixed(3)}`;
                }

                slider.addEventListener('pointerdown', () => { scrubbing = true; });
                const finishScrub = () => {
                    scrubbing = false;
                    previous = performance.now();
                };
                slider.addEventListener('pointerup', finishScrub);
                slider.addEventListener('pointercancel', finishScrub);
                slider.addEventListener('lostpointercapture', finishScrub);
                slider.addEventListener('input', () => {
                    phase = slider.valueAsNumber / 1000.0;
                    setTime(phase);
                });

                function animate(now) {
                    if (!scrubbing) {
                        phase += (now - previous) / (1000.0 * secondsPerRotation);
                        phase %= 1.0;
                        setTime(phase);
                    }
                    previous = now;
                    requestAnimationFrame(animate);
                }

                setTime(phase);
                requestAnimationFrame(animate);
            }""")
            return DOM.div(TIMESERIES_STYLE, root)
        end
    end
end

function export_regrid_timeseries(path = TIMESERIES_OUTPUT)
    mkpath(dirname(path))
    Bonito.export_static(path, regrid_timeseries_app())
    println("wrote ", path)
    return path
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_regrid_timeseries()
