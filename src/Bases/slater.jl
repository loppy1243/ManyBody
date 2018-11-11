export create, create!, annihil, annihil!

struct Slater{SPB<:TensorBasis} <: AbstractBasis
    bits::BitArray

    function Slater{SPB}(parts::BitArray) where B<:TensorBasis
        @assert size(parts) == fulldims(SPB)
        new(parts)
    end
    function Slater{SPB}(parts) where SPB<:TensorBasis
        bits = falses(fulldims(SPB))

        for p in parts
            p = convert(SPB, p)

            bits[index(p)] = true
        end

        new(bits)
    end
end
Slater{SPB}(ps...) where SPB<:TensorBasis = Slater{SPB}(ps)
Slater(ps::NTuple{<:Any, B}) where B<:TensorBasis = Slater{B}(ps)
Slater(ps::B...) where B<:TensorBasis = Slater{B}(ps)

spbasis(::Type{Slater{SPB}}) where SPB<:TensorBasis = SPB
spbasis(::Slater) = spbasis(typeof(Slater))

Base.:(==)(s1::B, s2::B) where B<:Slater = s1.bits == s2.bits

Base.in(p::TensorBasis, s::Slater) = s.bits[convert(spbasis(s), p)]

index(s::Slater) = sum(s.bits[i]*(1 << (i-1)) for i in LinearIndices(s.bits))
function indexbasis(B::Type{<:Slater}, ix::Int)
    SPB = spbasis(B)
    bits = falses(fulldims(SPB))
    for i in LinearIndices(bits)
        bits[i] = (ix & (1 << (i-1))) != 0
    end

    B(bits)
end

dim(B::Type{<:Slater}) = 2^dim(spbasis(B)) - 1
#nocc(s::Slater) = count(s.bits)
#nunocc(s::Slater) = count(.~s.bits)

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
