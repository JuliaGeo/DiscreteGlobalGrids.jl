# Denoise a Hikari frame with Intel Open Image Denoise.
#
# Hikari's own denoiser is an à-trous wavelet filter: a hand-tuned edge-stopping
# blur, and one that has to be kept to a single pass before it starts taking the
# terrain's ridges with the noise. OIDN is a trained network for exactly this
# job — a path-traced frame plus its first-hit albedo and normal — so it removes
# more noise per sample and keeps the detail the filter cannot distinguish from
# it.
#
# `Film` already carries everything the network wants: `framebuffer` is linear
# HDR sRGB, `albedo` and `normal` are the first-hit buffers the aux kernel
# writes. Nothing needs rendering twice.
#
# Install the release binaries and point `OIDN_ROOT` at them:
#
#     mkdir -p ~/.local/opt && cd ~/.local/opt
#     curl -sLO https://github.com/RenderKit/oidn/releases/download/v2.5.1/oidn-2.5.1.x86_64.linux.tar.gz
#     tar xzf oidn-2.5.1.x86_64.linux.tar.gz
#
# Then `include` this beside the scene and call `use_oidn!()`, which substitutes
# the network for `Hikari.denoise!`. Everything downstream — `--denoise`, the
# postprocess, `record_longrunning` — is unchanged.

module OIDNDenoise

using Hikari
using Colors: RGB
using GeometryBasics: Vec3f

const ROOT = get(ENV, "OIDN_ROOT",
    joinpath(homedir(), ".local", "opt", "oidn-2.5.1.x86_64.linux"))
const LIB = joinpath(ROOT, "lib", "libOpenImageDenoise.so")

const FORMAT_FLOAT3 = Cint(3)
const DEVICE_DEFAULT = Cint(0)   # CPU unless a GPU device library loads
const DEVICE_CPU = Cint(1)

# One device and one filter for the whole session. Building either costs far
# more than running it — the filter holds the network weights — and a video
# denoises the same shape every frame.
const DEVICE = Ref{Ptr{Cvoid}}(C_NULL)
const FILTER = Ref{Ptr{Cvoid}}(C_NULL)
const FILTER_SHAPE = Ref{Tuple{Int,Int}}((0, 0))

function device()
    if DEVICE[] == C_NULL
        isfile(LIB) || error("no OpenImageDenoise at $LIB; set OIDN_ROOT to the \
            unpacked release, or download it — see the top of this file")
        d = @ccall LIB.oidnNewDevice(DEVICE_DEFAULT::Cint)::Ptr{Cvoid}
        d == C_NULL && error("oidnNewDevice returned nothing")
        @ccall LIB.oidnCommitDevice(d::Ptr{Cvoid})::Cvoid
        checkerror(d)
        DEVICE[] = d
    end
    return DEVICE[]
end

function checkerror(d::Ptr{Cvoid})
    msg = Ref{Ptr{Cchar}}(C_NULL)
    code = @ccall LIB.oidnGetDeviceError(d::Ptr{Cvoid}, msg::Ptr{Ptr{Cchar}})::Cint
    code == 0 && return nothing
    text = msg[] == C_NULL ? "" : unsafe_string(msg[])
    error("OpenImageDenoise error $code: $text")
end

"""
    devicename() -> String

Which backend the network will run on. The release ships CPU, CUDA, HIP and SYCL
device libraries and loads whichever the machine can use, so this is worth
printing once rather than assuming.
"""
function devicename()
    d = device()
    n = @ccall LIB.oidnGetDeviceInt(d::Ptr{Cvoid}, "type"::Cstring)::Cint
    return get(Dict(1 => "CPU", 2 => "SYCL", 3 => "CUDA", 4 => "HIP", 5 => "Metal"),
               Int(n), "type $n")
end

# The network runs on whichever device OIDN picks, and on a GPU device it cannot
# reach host memory — `oidnSetSharedFilterImage` on a plain `Array` fails with
# "image data not accessible by the device". So every image goes through an
# `OIDNBuffer`, which the device allocates and which host writes and reads
# reach. The buffers are cached with the filter, since a video sends the same
# shape through every frame.
mutable struct Image
    buffer::Ptr{Cvoid}
    bytes::Int
end

const IMAGES = Dict{String,Image}()

function image(name::String, bytes::Int)
    img = get(IMAGES, name, nothing)
    if img === nothing || img.bytes != bytes
        img === nothing || @ccall LIB.oidnReleaseBuffer(img.buffer::Ptr{Cvoid})::Cvoid
        b = @ccall LIB.oidnNewBuffer(device()::Ptr{Cvoid}, bytes::Csize_t)::Ptr{Cvoid}
        b == C_NULL && error("oidnNewBuffer($bytes) returned nothing")
        img = Image(b, bytes)
        IMAGES[name] = img
    end
    return img
end

# A Julia matrix is column-major, so its columns are what OIDN should read as
# rows: passing `width = size(A, 1)` hands the network the image transposed,
# which costs nothing. Denoising is not orientation-sensitive and every buffer
# is transposed the same way. Normals are world-space vectors, so transposing
# the image does not touch them.
function setimage!(f::Ptr{Cvoid}, name::String, A::Matrix{T}) where {T}
    stride = sizeof(T)
    bytes = stride * length(A)
    img = image(name, bytes)
    @ccall LIB.oidnWriteBuffer(img.buffer::Ptr{Cvoid}, 0::Csize_t, bytes::Csize_t,
                               A::Ptr{Cvoid})::Cvoid
    @ccall LIB.oidnSetFilterImage(
        f::Ptr{Cvoid}, name::Cstring, img.buffer::Ptr{Cvoid}, FORMAT_FLOAT3::Cint,
        size(A, 1)::Csize_t, size(A, 2)::Csize_t,
        0::Csize_t, stride::Csize_t, (stride * size(A, 1))::Csize_t)::Cvoid
    return img
end

function readimage!(A::Matrix{T}, img::Image) where {T}
    @ccall LIB.oidnReadBuffer(img.buffer::Ptr{Cvoid}, 0::Csize_t, img.bytes::Csize_t,
                              A::Ptr{Cvoid})::Cvoid
    return A
end

function filter(shape::Tuple{Int,Int})
    if FILTER[] == C_NULL || FILTER_SHAPE[] != shape
        FILTER[] == C_NULL ||
            @ccall LIB.oidnReleaseFilter(FILTER[]::Ptr{Cvoid})::Cvoid
        f = @ccall LIB.oidnNewFilter(device()::Ptr{Cvoid}, "RT"::Cstring)::Ptr{Cvoid}
        f == C_NULL && error("oidnNewFilter(\"RT\") returned nothing")
        FILTER[] = f
        FILTER_SHAPE[] = shape
    end
    return FILTER[]
end

"""
    denoise!(color, albedo, normal) -> color

Run the network over `color` in place, guided by `albedo` and `normal`.

`hdr` is on because `framebuffer` is linear radiance with no tone map applied
yet, and `cleanAux` because the aux buffers come from a separate first-hit pass
rather than from the noisy path samples.
"""
function denoise!(color::Matrix{RGB{Float32}}, albedo::Matrix{RGB{Float32}},
                  normal::Matrix{Vec3f})
    size(color) == size(albedo) == size(normal) ||
        throw(DimensionMismatch("colour $(size(color)), albedo $(size(albedo)), \
            normal $(size(normal)) must agree"))
    f = filter(size(color))
    setimage!(f, "color", color)
    setimage!(f, "albedo", albedo)
    setimage!(f, "normal", normal)
    # Output shares the colour buffer: the network reads and writes the same
    # image, which OIDN supports and which saves an allocation the size of a
    # frame.
    out = image("color", sizeof(eltype(color)) * length(color))
    @ccall LIB.oidnSetFilterImage(
        f::Ptr{Cvoid}, "output"::Cstring, out.buffer::Ptr{Cvoid}, FORMAT_FLOAT3::Cint,
        size(color, 1)::Csize_t, size(color, 2)::Csize_t, 0::Csize_t,
        sizeof(eltype(color))::Csize_t,
        (sizeof(eltype(color)) * size(color, 1))::Csize_t)::Cvoid
    @ccall LIB.oidnSetFilterBool(f::Ptr{Cvoid}, "hdr"::Cstring, true::Bool)::Cvoid
    @ccall LIB.oidnSetFilterBool(f::Ptr{Cvoid}, "cleanAux"::Cstring, true::Bool)::Cvoid
    @ccall LIB.oidnCommitFilter(f::Ptr{Cvoid})::Cvoid
    @ccall LIB.oidnExecuteFilter(f::Ptr{Cvoid})::Cvoid
    checkerror(device())
    readimage!(color, out)
    return color
end

"""
    denoise_film!(film)

Denoise a `Hikari.Film` in place.

The buffers live on the GPU, so each one comes back to the host, through the
network, and returns. That round trip is the cost of using a library that does
not speak Vulkan: about 25 MB each way at 1080p.
"""
function denoise_film!(film::Hikari.Film)
    color = Array(film.framebuffer)
    albedo = Array(film.albedo)
    normal = Array(film.normal)
    denoise!(color, albedo, normal)
    copyto!(film.framebuffer, color)
    return film
end

"""
    use_oidn!()

Send `Hikari.denoise!` through the network instead of the à-trous filter.

`config` is accepted and ignored: the network has no iteration count or
edge-stopping widths to set, which is the point of it.
"""
function use_oidn!()
    @eval Hikari function denoise!(film::Hikari.Film; config = nothing)
        return $(@__MODULE__).denoise_film!(film)
    end
    @info "denoising with Open Image Denoise" device = devicename()
    return nothing
end

end # module
