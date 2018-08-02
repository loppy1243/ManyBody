export addpart!, addhole!, rmpart!, rmhole!

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

Base.in(p::SP, s::PartHole{R}) where {SP, R<:RefState{SP}} =
    ispart(R, p) ? s.parts[pindex(R, p)] : !s.holes[hindex(R, p)]

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

function addpart!(s::PartHole{R}, p::SP) where {SP, R<:RefState{SP}}
    i = pindex(R, p)
    ispart(R, p) && (s.parts[i] = true)

    1 - 2((count(.~s.holes) + count(s.parts[1:i])) % 2)
end
function addhole!(s::PartHole{R}, p::SP) where {SP, R<:RefState{SP}}
    i = hindex(R, p)
    ret = ishole(R, p) && (s.holes[i] = true)

    1 - 2((count(.~s.holes[1:i])) % 2)
end
function rmpart!(s::PartHole{R}, p::SP) where {SP, R<:RefState{SP}}
    i = pindex(R, p)
    ispart(R, p) && (s.parts[i] = false)

    1 - 2((count(.~s.holes) + count(s.parts[1:i])) % 2)
end
function rmhole!(s::PartHole{R}, p::SP) where {SP, R<:RefState{SP}}
    i = hindex(R, p)
    ishole(R, p) && (s.holes[i] = false)

    1 - 2((count(.~s.holes[1:i])) % 2)
end

RefStates.nholes(s::PartHole) = count(s.holes)
RefStates.nparts(s::PartHole) = count(s.parts)
