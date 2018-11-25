struct Slater{SPB<:AbstractBasis, N} <: MBBasis{SPB}
    bits::BitArray{N}

    function Slater{SPB, N}(parts::BitArray{N}) where {SPB<:TensorBasis, N}
        @assert size(parts) == fulldims(SPB)
        new(parts)
    end
    function Slater{SPB, 1}(parts::BitVector) where SPB<:AbstractBasis
        @assert length(parts) == dim(SPB)
        new(parts)
    end
end
_bits_size(SPB::Type{<:AbstractBasis}) = (dim(SPB),)
_bits_size(SPB::Type{<:TensorBasis}) = fulldims(SPB)
function Slater{SPB}(parts) where SPB<:AbstractBasis
    bits = falses(_bits_size(SPB))

    for p in parts
        p = convert(SPB, p)

        bits[p] = true
    end

    Slater{SPB, ndims(bits)}(bits)
    
end
Slater{SPB}(parts::BitArray) where SPB<:TensorBasis = Slater{SPB, ndims(parts)}(parts)
Slater{SPB}(ps...) where SPB<:AbstractBasis = Slater{SPB}(ps)
Slater(ps::NTuple{<:Any, SPB}) where SPB<:AbstractBasis = Slater{SPB}(ps)
Slater(ps::SPB...) where SPB<:AbstractBasis = Slater{SPB}(ps)
Slater(ref::RefState) = Slater{spbasis(ref)}(occ(ref))
Slater{SPB}(ref::RefState{SPB}) where SPB<:AbstractBasis = Slater{SPB}(occ(ref))

Base.copy(s::Slater) = typeof(s)(copy(s.bits))

### MBBasis methods
##############################################################################################
occ(s::Slater) = (SPI=indexer(spbasis(s)); map(i -> SPI[i], findall(s.bits)))
occinds(s::Slater) = findall(s.bits)
unocc(s::Slater) = (SPI=indexer(spbasis(s)); map(i -> SPI[i], findall(.~s.bits)))
unoccinds(s::Slater) = findall(.~s.bits)
nocc(s::Slater) = count(s.bits)
nunocc(s::Slater) = count(.~s.bits)
isocc(s::Slater, p) = s.bits[convert(spbasis(s), p)]

### AbstractBasis methods
##############################################################################################
dim(B::Type{<:Slater}) = 2^dim(spbasis(B)) - 1
index(s::Slater) = sum(s.bits[i]*(1 << (i-1)) for i in LinearIndices(s.bits))
function indexbasis(B::Type{<:Slater}, ix::Int)
    bits = falses(_bits_size(spbasis(B)))
    for i in LinearIndices(bits)
        bits[i] = (ix & (1 << (i-1))) != 0
    end

    B(bits)
end

### Slater methods
##############################################################################################
Base.:(==)(s1::B, s2::B) where B<:Slater = all(s1.bits .== s2.bits)

"""
    acsgn(s, p)

Return the sign picked up by applying an annihilation/creation operator of `p` to `s`.

Only returns +1 or -1.
"""
function acsgn(s::Slater{SPB}, p::SPB) where SPB<:AbstractBasis
    i = linearindex(p)
    1 - 2(count(s.bits[1:i-1]) % 2)
end

"""
    create!(s, p)

Apply a creation operator of `p` to `s`, changing `s` in place.

Returns `(sgn, occ)` where `sgn` is the acquired sign and `occ` is the original occupation of
`p` in `s`.
"""
function create!(s::Slater{SPB}, p::SPB) where SPB<:AbstractBasis
    sgn = acsgn(s, p)
    bit = s.bits[p]
    bit && return (0, bit)
    s.bits[p] = true
    (sgn, bit)
end


"""
    annihil!(s, p)

Apply a annihilation operator of `p` to `s`, changing `s` in place.

Returns `(sgn, occ)` where `sgn` is the acquired sign and `occ` is the original occupation of
`p` in `s`.
"""
function annihil!(s::Slater{SPB}, p::SPB) where SPB<:AbstractBasis
    p::spbasis(s) = p

    sgn = acsgn(s, p)
    bit = s.bits[p]
    ~bit && return (0, bit)
    s.bits[p] = false
    (sgn, bit)
end

function Base.show(io::IO, x::Slater)
    SPB = spbasis(x); SPI = indexer(SPB)
    str = prod((SPI[I] in x ? string(SPI[I]) : "_") * (I == lastindex(SPI) ? "" : ", ")
               for I in SPI)
    print(io, "Slater($(dim(SPB)))[$str]")
end

function Base.show(io::IO, ::MIME"text/plain", x::Slater{<:Pairing})
    SPB = spbasis(x)
    L = nlevels(SPB)
    for i = L:-1:1
        print(io, lpad(i, ndigits(L)), "|  ")
        if SPB(i, SPINDOWN) in x
            print(io, "↓  ")
        else
            print(io, "o  ")
        end

        if SPB(i, SPINUP) in x
            println(io, "↑")
        else
            println(io, "o")
        end
    end
end
