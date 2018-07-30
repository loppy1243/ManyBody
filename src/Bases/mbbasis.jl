struct PartHole{R<:RefState} <: MBBasis
    parts::BitVector
    holes::BitVector

    PartHole{R}(parts, holes) where R = new(parts, holes)
    function PartHole{R}(phs) where R
        ps = falses(nparts(R))
        hs = falses(nholes(R))

        for (p, h) in phs
            @assert ispart(R, p)
            @assert ishole(R, h)

            pn = pindex(R, p)
            hn = hindex(R, h)

            ps[pn] = true
            hs[hn] = true
        end

        new(ps, hs)
    end
end
PartHole{R}(phs::Vararg{Tuple{SP, SP}}) where {SP, R<:RefState{SP}} = PartHole{R}(phs)

Base.:(==)(s1::MB, s2::MB) where MB<:PartHole =
    s1.parts == s2.parts && s1.holes == s2.holes

function index(s::PartHole)
    bs = [s.holes; s.parts]

    sum(bs[i]*2^(i-1) for i in indices(bs, 1))
end

function indexbasis(::Type{PartHole{R}}, ix::Int) where R
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

    PartHole{R}(parts, holes)
end

dim(::Type{PartHole{R}}) where R = sum(0:n_occ(R)) do i
    binomial(n_occ(R), i)*binomial(n_unocc(R), i)
end
