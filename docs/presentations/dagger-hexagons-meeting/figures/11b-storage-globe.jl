isdefined(@__MODULE__, :DaggerTalkFigures) ||
    include(joinpath(@__DIR__, "11-dagger-theme.jl"))

using .DGGSTalkFigures, .DaggerTalkFigures
using Bonito, DiscreteGlobalGrids, DiscreteGlobalGridsVisualization
using FlyThroughPaths, LinearAlgebra, Makie, WGLMakie

const STORAGE_GLOBE_OUTPUT = joinpath(@__DIR__, "html", "11b-storage-globe.html")

unit_vector(lon, lat) = begin
    lambda, phi = deg2rad.((lon, lat))
    [cos(phi) * cos(lambda), cos(phi) * sin(lambda), sin(phi)]
end

local_north(lon, lat) = begin
    lambda, phi = deg2rad.((lon, lat))
    normalize([-sin(phi) * cos(lambda), -sin(phi) * sin(lambda), cos(phi)])
end

function storage_camera_path()
    far_direction = unit_vector(-80.0, 24.0)
    root_direction = unit_vector(-105.0, 40.0)
    up = local_north(-105.0, 40.0)

    far = ViewState(; eyeposition = 3.0far_direction,
        lookat = zeros(3), upvector = up, fov = 45.0)
    turned = ViewState(; eyeposition = 2.8root_direction,
        lookat = zeros(3), upvector = up, fov = 45.0)
    chunk = ViewState(; eyeposition = 1.055root_direction,
        lookat = root_direction, upvector = up, fov = 42.0)
    # Keep the selected parent fully in frame at maximum zoom.  Centering this
    # state on a leaf boundary caused the parent to slide off-canvas; the leaf
    # curves are already legible at this centered scale.
    detail = ViewState(; eyeposition = 1.045root_direction,
        lookat = root_direction, upvector = up, fov = 36.0)

    Path(far) * Pause(0.55) *
        ConstrainedMove(1.15, turned; constraint = :rotation, speed = :sinusoidal) *
        ConstrainedMove(2.0, chunk; constraint = :none, speed = :sinusoidal) *
        Pause(0.65) *
        ConstrainedMove(1.7, detail; constraint = :none, speed = :sinusoidal) *
        Pause(0.9)
end

function browser_camera_keyframes(pathspec; fps = 60)
    times = range(0, FlyThroughPaths.duration(pathspec);
        length = ceil(Int, FlyThroughPaths.duration(pathspec) * fps) + 1)
    [(eye = collect(state.eyeposition), target = collect(state.lookat),
        up = collect(state.upvector), fov = state.fov) for state in pathspec.(times)]
end

const STORAGE_GLOBE_STYLE = Styles(
    CSS("html, body", "width" => "100%", "height" => "100%", "margin" => "0",
        "display" => "grid", "place-items" => "center", "overflow" => "hidden",
        "background" => "#ffffff"),
    CSS(".storage-globe", "position" => "relative", "width" => "960px",
        "height" => "540px", "overflow" => "hidden", "background" => "#ffffff"),
    CSS(".storage-globe__plot", "position" => "absolute", "inset" => "0"),
    CSS(".storage-globe__timeline", "position" => "absolute", "z-index" => "4",
        "left" => "180px", "right" => "180px", "bottom" => "18px",
        "height" => "24px", "display" => "flex", "align-items" => "center",
        "padding" => "0 10px", "border" => "1px solid #dee2e6",
        "border-radius" => "2px", "background" => "rgba(255,255,255,0.90)"),
    CSS(".storage-globe__slider", "width" => "100%", "accent-color" => "#389826"),
)

function storage_globe_app()
    system = IGeo7System()
    chunks = levelgrid(system, 5)
    root = cellat(chunks, -105.0, 40.0)
    adjacent = collect(neighbors(chunks, root))
    boundary_cells = collect(border(subtree(system, root, 12); cells = true))
    adjacent_boundaries = [collect(border(subtree(system, neighbor, 12); cells = true))
        for neighbor in adjacent]

    pathspec = storage_camera_path()
    keyframes = browser_camera_keyframes(pathspec)

    with_theme(JG_THEME) do
        fig, body = slide_figure()
        axis = globe_axis(body[1, 1]; camera_longlat = (-80, 24),
            camera_altitude = 1.8)

        dggpoly!(axis, chunks; color = JG.paper_off,
            strokecolor = (JG.green_dark, 0.20), strokewidth = 0.22)
        coastlines!(axis)
        dggpoly!(axis, system, [root]; color = (JG.purple_100, 0.68),
            strokecolor = JG.purple, strokewidth = 1.8, zlevel = 0.021)
        for cells in adjacent_boundaries
            dggpoly!(axis, system, cells; color = (JG.paper, 0.0),
                strokecolor = (JG.green_dark, 0.94), strokewidth = 0.78,
                zlevel = 0.021)
        end
        dggpoly!(axis, system, boundary_cells; color = (JG.paper, 0.02),
            strokecolor = (JG.purple, 0.98), strokewidth = 0.86, zlevel = 0.021)

        slider = DOM.input(; class = "storage-globe__slider", type = "range",
            min = "0", max = "1000", step = "1", value = "0",
            ariaLabel = "Camera timeline")
        root_element = DOM.div(
            DOM.div(WGLMakie.WithConfig(fig; resize_to = :parent);
                class = "storage-globe__plot"),
            DOM.div(slider; class = "storage-globe__timeline");
            class = "storage-globe",
        )

        App() do session
            # Static Bonito exports have no Julia process.  WGLMakie's offline
            # OrbitControls therefore own the camera; a browser animation is
            # both the autoplay mechanism and the slider's scrub target.
            Bonito.onload(session, root_element, js"""async (root) => {
                const scene = await $(axis.scene);
                const controls = scene.orbitcontrols;
                const camera = controls && controls.object;
                const slider = root.querySelector('.storage-globe__slider');
                const keyframes = $(keyframes);

                if (!controls || !camera || !slider) {
                    console.error('storage globe: WGLMakie browser camera unavailable');
                    return;
                }

                let phase = 0.0;
                let previous = performance.now();
                let scrubbing = false;
                let running = true;
                const secondsPerRun = 8.0;
                const mix = (a, b, t) => a + (b - a) * t;
                const mix3 = (a, b, t) => [mix(a[0], b[0], t),
                    mix(a[1], b[1], t), mix(a[2], b[2], t)];

                function setCamera(u) {
                    const clamped = Math.min(1.0, Math.max(0.0, u));
                    const p = clamped * (keyframes.length - 1);
                    const i = Math.min(keyframes.length - 2, Math.floor(p));
                    const t = p - i;
                    const a = keyframes[i];
                    const b = keyframes[i + 1];
                    const eye = mix3(a.eye, b.eye, t);
                    const target = mix3(a.target, b.target, t);
                    const up = mix3(a.up, b.up, t);
                    camera.position.set(...eye);
                    camera.up.set(...up).normalize();
                    controls.target.set(...target);
                    camera.fov = mix(a.fov, b.fov, t);
                    camera.updateProjectionMatrix();
                    controls.update();
                    slider.value = String(Math.round(clamped * 1000));
                }

                slider.addEventListener('pointerdown', () => {
                    scrubbing = true;
                    running = false;
                });
                slider.addEventListener('pointerup', () => {
                    scrubbing = false;
                    previous = performance.now();
                    if (phase < 1.0 && !running) {
                        running = true;
                        requestAnimationFrame(animate);
                    }
                });
                slider.addEventListener('input', () => {
                    phase = slider.valueAsNumber / 1000.0;
                    setCamera(phase);
                });

                function animate(now) {
                    if (!scrubbing) {
                        phase = Math.min(1.0, phase + (now - previous) /
                            (1000.0 * secondsPerRun));
                        setCamera(phase);
                    }
                    previous = now;
                    if (phase < 1.0 && !scrubbing) {
                        requestAnimationFrame(animate);
                    } else {
                        running = false;
                    }
                }

                setCamera(phase);
                requestAnimationFrame(animate);
            }""")
            DOM.div(STORAGE_GLOBE_STYLE, root_element)
        end
    end
end

function export_storage_globe(path = STORAGE_GLOBE_OUTPUT)
    mkpath(dirname(path))
    Bonito.export_static(path, storage_globe_app())
    println("wrote ", path)
    return path
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_storage_globe()
