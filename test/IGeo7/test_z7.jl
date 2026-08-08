# Tests for the pure-integer Z7 layer (`src/IGeo7/z7.jl`).
#
# Authorities:
#   - spec/design.md §3 (cell-id scheme) and §4.3 (hierarchy ops)
#   - spec/interface-contract.md ("Z7 index", "Native API")
#   - spec/aperture7-indexing-spec.md §4 (bit format, prefix ops, counting)
#   - test/IGeo7/vectors/pentagon_chains.csv, test/IGeo7/vectors/hierarchy.csv

using Test

# `_Z7` gives access to the names `IGeo7` does not export (`z7_deleted_digit`,
# `Z7_MAX_RESOLUTION`, ...); the exported ones come from `runtests.jl`'s
# `using .IGeo7`.
const _Z7 = IGeo7

const _Z7_VECTORS = joinpath(@__DIR__, "vectors")

"""Measure allocations of `f(args...)` after a warm-up call, behind a function
barrier so the measurement is not polluted by boxing of untyped globals."""
function _z7_alloc(f::F, args...) where {F}
    f(args...)
    return @allocated f(args...)
end

"""Read a headered CSV into a header vector and a vector of row vectors."""
function _z7_read_csv(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    rows = [split(line, ',') for line in lines[2:end] if !isempty(line)]
    return header, rows
end

@testset "Z7 bit layer" begin

    @testset "constants and layout" begin
        @test MAX_RESOLUTION == 19
        @test _Z7.Z7_MAX_RESOLUTION == 20
        # Deleted-digit table [contract]: bases 0..5 drop 2, bases 6..11 drop 5.
        @test [_Z7.z7_deleted_digit(b) for b in 0:11] == [2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5]
        @test_throws InvalidZ7Error _Z7.z7_deleted_digit(12)
        @test_throws InvalidZ7Error _Z7.z7_deleted_digit(-1)
    end

    @testset "published bit-format example" begin
        # [contract] z7_from_string("0800433") == 0x80237fffffffffff
        published = z7_from_string("0800433")
        @test published === 0x80237fffffffffff
        @test z7_to_string(published) == "0800433"
        @test z7_to_hex(published) == "80237fffffffffff"
        @test z7_from_hex("80237fffffffffff") === published
        @test z7_from_hex("0x80237fffffffffff") === published
        @test z7_from_hex("0X80237FFFFFFFFFFF") === published
        @test z7_base_cell(published) == 8
        @test z7_resolution(published) == 5
        @test z7_digit(published, 1) == 0
        @test z7_digit(published, 2) == 0
        @test z7_digit(published, 3) == 4
        @test z7_digit(published, 4) == 3
        @test z7_digit(published, 5) == 3
        @test z7_digit(published, 6) == 7
        @test z7_digit(published, 20) == 7
        @test !z7_is_pentagon(published)
        @test is_valid_z7(published)

        parent = z7_parent(published)
        @test z7_to_string(parent) == "080043"
        @test z7_parent(published, 2) === z7_from_string("0800")
        @test z7_parent(published, 5) === published
        @test z7_is_descendant(published, parent)
        @test !z7_is_descendant(parent, published)
    end

    @testset "base-cell strings and res-0 pentagons" begin
        @test z7_from_string("00") === 0x0fffffffffffffff
        @test z7_from_string("11") === 0xbfffffffffffffff
        for base in 0:11
            text = lpad(string(base), 2, '0')
            z = z7_from_string(text)
            @test z === (UInt64(base) << 60) | 0x0fffffffffffffff
            @test z7_base_cell(z) == base
            @test z7_resolution(z) == 0
            @test z7_is_pentagon(z)
            @test is_valid_z7(z)
            @test z7_to_string(z) == text
            @test z7_from_hex(z7_to_hex(z)) === z
        end
        # ascending base order == ascending u64
        ids = [z7_from_string(lpad(string(b), 2, '0')) for b in 0:11]
        @test issorted(ids) && allunique(ids)
    end

    @testset "children of pentagons and hexagons" begin
        north = z7_from_string("00")
        @test length(z7_children(north)) == 6
        @test z7_to_string.(collect(z7_children(north))) ==
              ["000", "001", "003", "004", "005", "006"]
        @test_throws InvalidZ7Error z7_child(north, 2)

        south = z7_from_string("06")
        @test length(z7_children(south)) == 6
        @test z7_to_string.(collect(z7_children(south))) ==
              ["060", "061", "062", "063", "064", "066"]
        @test_throws InvalidZ7Error z7_child(south, 5)

        hexagon = z7_from_string("001")
        @test length(z7_children(hexagon)) == 7
        @test z7_to_string.(collect(z7_children(hexagon))) ==
              ["0010", "0011", "0012", "0013", "0014", "0015", "0016"]
        # A hexagon prefix keeps every digit, including the base's deleted one.
        @test z7_to_string(z7_child(hexagon, 2)) == "0012"

        # Center child of a pentagon is a pentagon; the others are hexagons.
        @test z7_is_pentagon(z7_child(north, 0))
        @test !z7_is_pentagon(z7_child(north, 1))

        # Children are ascending and their parent is the original id.
        for z in (north, south, hexagon)
            kids = collect(z7_children(z))
            @test issorted(kids)
            @test all(z7_parent(k) === z for k in kids)
            @test all(z7_is_descendant(k, z) for k in kids)
            @test all(z7_resolution(k) == z7_resolution(z) + 1 for k in kids)
        end

        @test_throws InvalidZ7Error z7_child(north, 7)
        @test_throws InvalidZ7Error z7_child(north, -1)
        @test_throws InvalidZ7Error z7_child(north, 8)
    end

    @testset "string parsing and validation" begin
        @test_throws InvalidZ7Error z7_from_string("")
        @test_throws InvalidZ7Error z7_from_string("0")
        @test_throws InvalidZ7Error z7_from_string("12")
        @test_throws InvalidZ7Error z7_from_string("99")
        @test_throws InvalidZ7Error z7_from_string("002")   # deleted digit, base 00
        @test_throws InvalidZ7Error z7_from_string("065")   # deleted digit, base 06
        @test_throws InvalidZ7Error z7_from_string("0002")  # deleted deeper in the chain
        @test_throws InvalidZ7Error z7_from_string("0107")  # digit 7 is padding only
        @test_throws InvalidZ7Error z7_from_string("0a")
        @test_throws InvalidZ7Error z7_from_string("0x00")
        @test_throws InvalidZ7Error z7_from_string("01" * repeat("0", 21))  # 21 digits

        # Deleted digit only applies while the prefix is all-zero.
        @test z7_to_string(z7_from_string("0012")) == "0012"
        @test z7_to_string(z7_from_string("00102")) == "00102"
        @test_throws InvalidZ7Error z7_from_string("00002")
        @test_throws InvalidZ7Error z7_from_string("000002")

        # Malformed padding: a non-7 slot after the first 7 slot.
        malformed = z7_from_string("00") & ~(UInt64(7) << 54)
        @test !is_valid_z7(malformed)
        @test_throws InvalidZ7Error z7_resolution(malformed)
        @test_throws InvalidZ7Error z7_to_string(malformed)

        # Out-of-range base cell in the raw bits.
        @test !is_valid_z7(0xcfffffffffffffff)
        @test !is_valid_z7(0xffffffffffffffff)

        # Deleted-pentagon subsequence in the raw bits ("001" -> "002").
        @test !is_valid_z7(z7_from_string("001") ⊻ (UInt64(3) << 57))
    end

    @testset "hex parsing" begin
        @test_throws InvalidZ7Error z7_from_hex("80237fffffffffx")
        @test_throws InvalidZ7Error z7_from_hex("80237fffffffff")     # 14 chars
        @test_throws InvalidZ7Error z7_from_hex("080237fffffffffff")  # 17 chars
        @test_throws InvalidZ7Error z7_from_hex("0x")
        @test_throws InvalidZ7Error z7_from_hex(" 80237ffffffffff")   # padded, 15 digits
        @test_throws InvalidZ7Error z7_from_hex("cfffffffffffffff")   # base 12
        @test_throws InvalidZ7Error z7_from_hex("0e3fffffffffffff")   # bad padding
        @test z7_from_string("00") & ~(UInt64(7) << 54) === 0x0e3fffffffffffff
        @test all(c -> c in "0123456789abcdef", z7_to_hex(z7_from_string("0800433")))
        @test length(z7_to_hex(z7_from_string("00"))) == 16
    end

    @testset "resolution 20 prefix ops" begin
        text20 = "01" * repeat("0", 20)
        z20 = z7_from_string(text20)
        @test z7_resolution(z20) == 20
        @test z7_to_string(z20) == text20
        @test is_valid_z7(z20)
        @test z7_is_pentagon(z20)
        @test z7_from_hex(z7_to_hex(z20)) === z20
        @test z7_digit(z20, 20) == 0
        @test_throws InvalidZ7Error z7_digit(z20, 21)
        @test_throws InvalidZ7Error z7_digit(z20, 0)

        # Prefix ops stay valid at resolution 20 ...
        @test z7_resolution(z7_parent(z20)) == 19
        @test z7_to_string(z7_parent(z20)) == "01" * repeat("0", 19)
        @test z7_parent(z20, 0) === z7_from_string("01")
        @test z7_is_descendant(z20, z7_from_string("01"))
        # ... but there is no resolution 21.
        @test_throws InvalidZ7Error z7_children(z20)
        @test_throws InvalidZ7Error z7_child(z20, 0)

        # A resolution-19 id still has children (they are resolution 20).
        z19 = z7_parent(z20)
        kids = z7_children(z19)
        @test length(kids) == 6                     # pentagon chain
        @test all(z7_resolution(k) == 20 for k in kids)
        @test z7_child(z19, 1) === kids[2]

        # Deepest hexagon: 7 children at resolution 20.
        deep = z7_from_string("01" * "1" * repeat("0", 18))
        @test z7_resolution(deep) == 19
        @test length(z7_children(deep)) == 7
    end

    @testset "parent / descendant edge cases" begin
        root = z7_from_string("00")
        @test_throws InvalidZ7Error z7_parent(root)          # no parent at res 0
        @test_throws InvalidZ7Error z7_parent(root, 1)       # finer than the id
        @test_throws InvalidZ7Error z7_parent(root, -1)
        @test_throws InvalidZ7Error z7_parent(root, 21)
        @test z7_parent(root, 0) === root

        a = z7_from_string("0800433")
        @test z7_is_descendant(a, a)                        # reflexive [design §3]
        @test !z7_is_descendant(a, z7_from_string("0100433"))
        @test !z7_is_descendant(z7_from_string("080044"), z7_from_string("080043"))
        for r in 0:5
            @test z7_is_descendant(a, z7_parent(a, r))
        end
    end

    @testset "digit accessors" begin
        z = z7_from_string("0054321")
        @test [z7_digit(z, k) for k in 1:7] == [5, 4, 3, 2, 1, 7, 7]
        @test z7_base_cell(z) == 0
        @test z7_resolution(z) == 5
        # digits round-trip through z7_child appends
        built = z7_from_string("00")
        for d in (5, 4, 3, 2, 1)
            built = z7_child(built, d)
        end
        @test built === z
        @test [z7_digit(built, k) for k in 1:5] == [5, 4, 3, 2, 1]
    end

    @testset "oracle vectors: pentagon_chains.csv" begin
        header, rows = _z7_read_csv(joinpath(_Z7_VECTORS, "pentagon_chains.csv"))
        @test header == ["base", "depth", "pentagon_z7_string", "children_z7_strings",
                         "missing_digits"]
        @test !isempty(rows)
        for row in rows
            base = parse(Int, row[1])
            depth = parse(Int, row[2])
            text = String(row[3])
            expected_children = String.(split(row[4], ';'))
            missing_digit = parse(Int, row[5])

            z = z7_from_string(text)
            @test z7_base_cell(z) == base
            @test z7_resolution(z) == depth - 1
            @test z7_is_pentagon(z)
            @test is_valid_z7(z)
            @test z7_to_string(z) == text
            @test _Z7.z7_deleted_digit(base) == missing_digit

            kids = z7_children(z)
            @test length(kids) == 6
            @test z7_to_string.(collect(kids)) == expected_children
            @test issorted(collect(kids))
            @test all(z7_parent(k) === z for k in kids)

            # The missing digit is rejected by both the child op and the parser.
            @test_throws InvalidZ7Error z7_child(z, missing_digit)
            @test_throws InvalidZ7Error z7_from_string(text * string(missing_digit))
            # Every other digit is accepted.
            for d in 0:6
                d == missing_digit && continue
                @test z7_to_string(z7_child(z, d)) == text * string(d)
            end
        end
    end

    @testset "oracle vectors: hierarchy.csv" begin
        header, rows = _z7_read_csv(joinpath(_Z7_VECTORS, "hierarchy.csv"))
        @test header == ["z7_hex", "z7_string", "parent_chain_strings", "children_strings"]
        @test !isempty(rows)
        for row in rows
            hex = String(row[1])
            text = String(row[2])
            parents = String.(split(row[3], ';'))
            expected_children = String.(split(row[4], ';'))

            z = z7_from_hex(hex)
            @test z === z7_from_string(text)
            @test z7_to_hex(z) == hex
            @test z7_to_string(z) == text
            @test is_valid_z7(z)

            res = length(text) - 2
            @test z7_resolution(z) == res
            @test z7_base_cell(z) == parse(Int, text[1:2])

            # Parent chain: coarser by one resolution per entry, res-1 down to 0.
            @test length(parents) == res
            @test [z7_to_string(z7_parent(z, r)) for r in (res-1):-1:0] == parents
            walker = z
            for expected in parents
                walker = z7_parent(walker)
                @test z7_to_string(walker) == expected
            end
            @test all(z7_is_descendant(z, z7_from_string(p)) for p in parents)
            @test all(!z7_is_descendant(z7_from_string(p), z) for p in parents)

            kids = collect(z7_children(z))
            @test z7_to_string.(kids) == expected_children
            @test issorted(kids)
            @test all(z7_parent(k) === z for k in kids)
            @test all(z7_is_descendant(k, z) for k in kids)
        end
    end

    @testset "subtree enumeration and counts" begin
        # Descendant counts [a7 §4.3]: 7^d below a hexagon, (5*7^d + 1)/6 below
        # a pentagon (its center child stays a pentagon, one digit is deleted).
        for base in 0:11
            root = z7_from_string(lpad(string(base), 2, '0'))
            level = [root]
            for d in 1:3
                level = reduce(vcat, (collect(z7_children(z)) for z in level))
                @test length(level) == (5 * 7^d + 1) ÷ 6
                @test issorted(level) && allunique(level)
                @test all(is_valid_z7, level)
                @test all(z -> z7_resolution(z) == d, level)
                @test all(z -> z7_is_descendant(z, root), level)
                @test count(z7_is_pentagon, level) == 1
                # enumeration order == digit-lexicographic string order
                @test issorted(z7_to_string.(level))
            end
        end

        hexagon = z7_from_string("001")
        level = [hexagon]
        for d in 1:3
            level = reduce(vcat, (collect(z7_children(z)) for z in level))
            @test length(level) == 7^d
            @test issorted(level) && allunique(level)
            @test all(z -> z7_is_descendant(z, hexagon), level)
            @test count(z7_is_pentagon, level) == 0
        end

        # World totals [contract]: num_cells(r) == 12 * p(r) == 10*7^r + 2.
        for r in 0:3
            @test 12 * ((5 * 7^r + 1) ÷ 6) == 10 * 7^r + 2
        end
    end

    @testset "type stability" begin
        z = z7_from_string("0800433")
        @test @inferred(z7_base_cell(z)) == 8
        @test @inferred(z7_resolution(z)) == 5
        @test @inferred(z7_digit(z, 3)) == 4
        @test @inferred(z7_is_pentagon(z)) === false
        @test @inferred(is_valid_z7(z)) === true
        @test @inferred(z7_parent(z)) === z7_from_string("080043")
        @test @inferred(z7_parent(z, 2)) === z7_from_string("0800")
        @test @inferred(z7_child(z, 1)) === z7_from_string("08004331")
        @test @inferred(z7_is_descendant(z, z7_parent(z))) === true
        @test length(@inferred(z7_children(z))) == 7
        @test isbitstype(typeof(z7_children(z)))
        @test @inferred(z7_to_string(z)) == "0800433"
        @test @inferred(z7_from_string("0800433")) === z
        @test @inferred(z7_to_hex(z)) == "80237fffffffffff"
        @test @inferred(z7_from_hex("80237fffffffffff")) === z
    end

    @testset "scalar paths are allocation-free" begin
        z = z7_from_string("0" * "8" * "00433" * repeat("1", 12))  # deep sample
        @test z7_resolution(z) == 17
        @test _z7_alloc(z7_resolution, z) == 0
        @test _z7_alloc(z7_base_cell, z) == 0
        @test _z7_alloc(z7_parent, z) == 0
        @test _z7_alloc(z7_parent, z, 3) == 0
        @test _z7_alloc(z7_child, z, 4) == 0
        @test _z7_alloc(z7_digit, z, 3) == 0
        @test _z7_alloc(z7_is_pentagon, z) == 0
        @test _z7_alloc(z7_is_descendant, z, z7_parent(z)) == 0
        @test _z7_alloc(is_valid_z7, z) == 0
        @test _z7_alloc(z7_children, z) == 0
        # the codecs validate without interpolating: parsing a *valid* input
        # allocates nothing even though every failure edge quotes the input
        @test _z7_alloc(z7_from_string, "0800433") == 0
        @test _z7_alloc(z7_from_hex, "80237fffffffffff") == 0
    end

    @testset "InvalidZ7Error: lazy messages" begin
        # The validation error carries only the facts; the sentence is built
        # in `showerror`, i.e. when the error is *printed*, so every throw
        # site is a plain `throw(InvalidZ7Error(...))` that stores immediates
        # and no formatted string (see src/z7.jl).
        @test InvalidZ7Error <: Exception
        @test fieldnames(InvalidZ7Error) == (:reason, :value, :got, :limit, :input)
        @test fieldtypes(InvalidZ7Error) == (Symbol, UInt64, Int, Int, String)
        # The type is concrete (no abstract/parametric fields), so a throw site
        # is a fixed-layout struct construction.
        @test isconcretetype(InvalidZ7Error)
        # Tradeoff: `input` is a `String` (so the offending text can be quoted
        # back) and `reason` is a `Symbol`, so the struct is not `isbits` — the
        # numeric payload is. That costs nothing on the non-throwing path (both
        # reference fields are set only on the cold edge that is about to throw),
        # which is what the allocation testset above asserts.
        @test !isbitstype(InvalidZ7Error)
        @test isbitstype(Tuple{UInt64,Int,Int})
        @test InvalidZ7Error(:base_range, zero(UInt64), 12, 11).input == ""

        # The premise of the design: constructing the error is much cheaper
        # than rendering its message.
        mk() = InvalidZ7Error(:base_range, zero(UInt64), 12, 11)
        e0 = mk()
        render(e) = sprint(showerror, e)
        render(e0)
        a_construct = @allocated mk()
        a_render = @allocated render(e0)
        @test a_construct < a_render

        root = z7_from_string("00")
        res20 = z7_from_string("01" * repeat("0", 20))
        malformed = z7_from_string("00") & ~(UInt64(7) << 54)

        # (thunk, reason, substrings the printed message must contain)
        thrown_cases = [
            (() -> _Z7.z7_deleted_digit(12), :base_range,
                ["base cell", "0:11", "12"]),
            (() -> z7_from_string("12"), :base_range,
                ["\"12\"", "base cell", "0:11"]),
            (() -> z7_resolution(malformed), :bad_padding,
                ["malformed", "0x0e3fffffffffffff", "padding"]),
            (() -> z7_from_hex("cfffffffffffffff"), :invalid_index,
                ["invalid Z7 index", "0xcfffffffffffffff", "pentagon"]),
            (() -> z7_child(root, 2), :deleted_digit,
                ["deleted", "pentagon chain", "0x0fffffffffffffff"]),
            (() -> z7_from_string("002"), :deleted_digit,
                ["\"002\"", "deleted", "pentagon chain"]),
            (() -> z7_digit(root, 21), :digit_level,
                ["digit level", "1:20", "21"]),
            (() -> z7_child(root, 7), :child_digit,
                ["child digit", "0:6", "7"]),
            (() -> z7_children(res20), :no_children,
                ["0x1000000000000000", "20", "children"]),
            (() -> z7_parent(root), :no_parent,
                ["0x0fffffffffffffff", "resolution-0", "parent"]),
            (() -> z7_parent(root, 3), :parent_res,
                ["ancestor resolution", "0:0", "3", "0x0fffffffffffffff"]),
            (() -> z7_from_string("0"), :bad_length,
                ["\"0\"", "0:20", "1 character"]),
            (() -> z7_from_string("0a"), :bad_base_char,
                ["\"0a\"", "'a'", "decimal"]),
            (() -> z7_from_string("0107"), :bad_digit_char,
                ["\"0107\"", "'7'", "0:6"]),
            (() -> z7_from_hex("0x"), :bad_hex,
                ["\"0x\"", "16 hexadecimal"]),
        ]
        for (thunk, reason, needles) in thrown_cases
            err = try
                thunk()
                nothing
            catch caught
                caught
            end
            @test err isa InvalidZ7Error
            @test err.reason === reason
            message = sprint(showerror, err)
            @test !isempty(message)
            @test !occursin('$', message)          # no botched interpolation
            @test all(occursin(n, message) for n in needles)
        end

        # The four grid-layer reasons belong to `src/IGeo7/grid.jl`, not to this
        # file's layer, so their messages are checked by constructing the error
        # directly; grid.jl's own suites cover the throw sites.
        built_cases = [
            (InvalidZ7Error(:res20_geometry, res20, 20, MAX_RESOLUTION),
                ["0x1000000000000000", "no geometry", "19"]),
            (InvalidZ7Error(:resolution_range, zero(UInt64), 20, MAX_RESOLUTION),
                ["resolution", "0:19", "20"]),
            (InvalidZ7Error(:descendant_res, root, 0, 1),
                ["descendant resolution", "1:19", "0x0fffffffffffffff"]),
            (InvalidZ7Error(:no_child_geometry, root, 0, MAX_RESOLUTION),
                ["0x0fffffffffffffff", "19", "no geometry"]),
        ]
        for (err, needles) in built_cases
            message = sprint(showerror, err)
            @test !occursin('$', message)
            @test all(occursin(n, message) for n in needles)
        end

        # Documented taxonomy, and nothing outside it in use.
        taxonomy = Set([:base_range, :bad_padding, :invalid_index, :deleted_digit,
            :digit_level, :child_digit, :no_children, :no_parent, :parent_res,
            :bad_length, :bad_base_char, :bad_digit_char, :bad_hex,
            :res20_geometry, :resolution_range, :descendant_res,
            :no_child_geometry])
        exercised = Set(vcat([c[2] for c in thrown_cases],
            [e.reason for (e, _) in built_cases]))
        @test exercised == taxonomy

        # An unrecognized reason still prints every field rather than failing.
        fallback = sprint(showerror, InvalidZ7Error(:nonesuch, UInt64(0x1234), 5, 6, "zz"))
        @test occursin("nonesuch", fallback)
        @test occursin("0x0000000000001234", fallback)
        @test occursin("zz", fallback)
    end
end
