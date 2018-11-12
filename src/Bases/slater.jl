struct Slater{SPB<:AbstractBasis} <: MBBasis{SPB}
    bits::BitArray

    function Slater{SPB}(parts::BitArray) where SPB<:TensorBasis
        @assert size(parts) == fulldims(SPB)
        new(parts)
    end
    function Slater{SPB}(parts::BitVector) where SPB<:AbstractBasis
        @assert length(parts) == dim(SPB)
        new(parts)
    end
end
@generated function Slater{SPB}(parts) where SPB<:AbstractBasis
    bits_size_expr = SPB<:TensorBasis ? :(fulldims(SPB)) : :(dim(SPB))

    quote
        bits = falses($bits_size_expr)

        for p in parts
            p = convert(SPB, p)

            bits[p] = true
        end

        Slater{SPB}(bits)
    end
end
Slater{SPB}(ps...) where SPB<:AbstractBasis = Slater{SPB}(ps)
Slater(ps::NTuple{<:Any, SPB}) where SPB<:AbstractBasis = Slater{SPB}(ps)
Slater(ps::SPB...) where SPB<:AbstractBasis = Slater{SPB}(ps)
Slater(ref::RefState) = Slater{spbasis(ref)}(occ(ref))
Slater{SPB}(ref::RefState{SPB}) where SPB<:AbstractBasis = Slater{SPB}(occ(ref))

### MBBasis methods
##############################################################################################
occ(s::Slater) = (SPB=spbasis(s); map(i -> SPB[i], findall(s.bits)))
unocc(s::Slater) = (SPB=spbasis(s); map(i -> SPB[i], findall(.~s.bits)))
nocc(s::Slater) = count(s.bits)
nunocc(s::Slater) = count(.~s.bits)
isocc(s::Slater, p) = s.bits[convert(spbasis(s), p)]

### AbstractBasis methods
##############################################################################################
dim(B::Type{<:Slater}) = 2^dim(spbasis(B)) - 1
index(s::Slater) = sum(s.bits[i]*(1 << (i-1)) for i in LinearIndices(s.bits))
@generated function indexbasis(::Type{B}, ix::Int) where B<:Slater
    bits_size_expr = spbasis(B)<:TensorBasis ? :(fulldims(SPB)) : :(dim(SPB))

    quote
        bits = falses($bits_size_expr)
        for i in LinearIndices(bits)
            bits[i] = (ix & (1 << (i-1))) != 0
        end

        B(bits)
    end
end

### Slater methods
##############################################################################################
Base.:(==)(s1::B, s2::B) where B<:Slater = s1.bits == s2.bits

function create!(s::Slater, p)
    p::spbasis(s) = p

    i = index(p)
    s.bits[i] && return 0
    s.bits[i] = true

    1 - 2(count(s.bits[1:i-1]) % 2)
end
create(s, p) = (s2 = deepcopy(s); (create!(s2, p), s2))

function annihil!(s::Slater, p)
    p::spbasis(s) = p

    i = index(p)
    ~s.bits[i] && return 0
    s.bits[i] = false

    1 - 2(count(s.bits[1:i-1]) % 2)
end
annihil(s, p) = (s2 = deepcopy(s); (annihil!(s2, p), s2))

function Base.show(io::IO, x::Slater)
    SPB = spbasis(x)
    str = prod((SPB[I] in x ? string(SPB[I]) : "_") * (I == lastindex(SB) ? "" : ", ")
               for I in eachindex(SPB))
    print(io, "Slater($(dim(SPB)))[$str]")
end
