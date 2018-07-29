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
