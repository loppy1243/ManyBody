struct MBBasis{R<:RefState, SP<:SPBasis} <: Basis
    parts::BitVector
    holes::BitVector

    MBBasis{R, SP}(parts, holes) where {R, SP} = new(parts, holes)
    function MBBasis{R, SP}(phs) where {R, SP}
        ps = falses(nholes(R))
        hs = falses(nparts(R))

        for (p, h) in phs
            pn = pindex(R, p)
            hn = hindex(R, h)
            @assert pn != hn
            @assert !ps[pn] && !hs[hn]

            ps[pn] = true
            hs[hn] = true
        end

        new(ps, hs)
    end
end
MBBasis{R, SP}(phs::Vararg{Tuple{SP, SP}}) where {R, SP} = MBBasis{R, SP}(phs)
MBBasis{R, SP}(phs::Vararg{Tuple{Int, Int}}) where {R, SP} = MBBasis{R, SP}(phs)
MBBasis(phs::Vararg{Tuple{SP, SP}}) where SP = MBBasis{Vaccuum{SP}, SP}(phs)

Base.:(==)(::MBBasis, ::MBBasis) = false
Base.:(==)(s1::MB, s2::MB) where MB<:MBBasis =
    s1.parts == s2.parts && s1.holes == s2.holes

function index(s::MBBasis)
    bs = [s.holes; s.parts]

    sum(bs[i]*2^(i-1) for i in indices(bs, 1))
end

function indexbasis(::Type{MBBasis{R, SP}}, ix::Int) where {R, SP}
    numhs = nholes(R)
    numps = nparts(R)

    holes = BitVector(numhs)
    parts = BitVector(numps)

    for i in 1:numhs
        holes[i] = (ix & (1 << (i-1))) != 0
    end
    for i in 1:numps
        parts[i] = (ix & (1 << (numhs + i - 1))) != 0
    end

    MBBasis{R, SP}(parts, holes)
end

dim(::Type{MBBasis{RS, SP}}) where {RS, SP} = 2^dim(SP) - 1


abstract type SubBasis
macro def_SubBasis(expr::Expr)
end

@def_SubBasis Paired{R, L} <: MBBasis{R, Pairing{L}}
struct PairedIndices{R, L, N}
    indices::NTuple{N, Int}

    PairedIndices{R, L, N}() where {R, L, N} = new([ixs(R, L, N)...])
end
PairedIndices{R, L}() where {R, L} = PairedIndices{R, L, x}()

struct Paired{R, L}
    index::Int

    Paired{R, L}(x::Int) = new(x)
    Paired{R, L}(x::MBBasis{R, Pairing{L}}) = new(index(x))
end

struct SubBasis{B<:Basis, N}
    indices::Vector{Int}

    function SubBasis{B, N}(ixs) where {B<:Basis, N}
        ixs = collect(Int, ixs) |> sort

        @assert ixs == unique(ixs)
        @assert length(ixs) == N
        @assert ixs[end] <= dim(B)

        new(ixs)
    end
end
SubBasis{B}(xs...) where B = SubBasis{B}(xs)
SubBasis{B}(xs::NTuple{N, B}) where {N, B<:Basis} = SubBasis{B, N}(map(index, xs))
SubBasis{B}(xs::NTuple{N, Int}) where {N, B} = SubBasis{B, N}(xs)

SubBasis(ixs...) = SubBasis(ixs)
SubBasis(ss::NTuple{N, B}) where {B<:Basis, N} = SubBasis{B, N}(map(index, ss))

dim(::Type{SubBasis{B, N}}) where {B, N} = N
