export create, create!, annihil, annihil!

struct Slater{B<:ConcreteBasis} <: ConcreteBasis
    bits::BitArray

    function Slater{B}(parts::BitArray) where B<:ConcreteBasis
        @assert size(parts) == innerdims(B)
        new(parts)
    end
    function Slater{B}(parts) where B<:ConcreteBasis
        bits = falses(innerdims(B))

        for p in parts
            p = convert(B, p)

            bits[index(p)] = true
        end

        new(bits)
    end
end
Slater{B}(ps...) where B<:ConcreteBasis = Slater{B}(ps)
Slater(ps::NTuple{<:Any, B}) where B<:ConcreteBasis = Slater{B}(ps)
Slater(ps::B...) where B<:ConcreteBasis = Slater{B}(ps)

IndexType(::Type{<:Slater}) = IndexTypes.Linear()

Base.:(==)(s1::B, s2::B) where B<:Slater = s1.bits == s2.bits

Base.in(p::MaybeIndex{B}, s::Slater{B}) where B<:ConcreteBasis = s.bits[index(p)]
Base.in(p::Rep, s::Rep) = inner(p) in inner(s)

innertype(::Type{Slater{B}}) where B<:ConcreteBasis = B
function inner(s::Slater)
    SB = innertype(s)
    SBI = indextype(B)
    Product(SB[I] for I in eachindex(SB) if SBI[I] in s)
end

index(s::Slater) = sum(s.bits[i]*2^(i-1) for i in LinearIndices(s.bits))
function indexbasis(B::Type{<:Slater}, ix::Int)
    SB = innertype(B)
    bits = falses(innerdims(SB))
    for i in LinearIndices(bits)
        bits[i] = (ix & (1 << (i-1))) != 0
    end

    Slater{B}(bits)
end

dim(::Type{Slater{B}}) where B = 2^dim(B) - 1
n_occ(s::Slater) = count(s.bits)
n_unoc(s::Slater) = count(.~s.bits)

create!(s::Slater{B}, p::Wrapped) where B = create!(s, inner(p))
function create!(s::Slater{B}, p::MaybeIndex{B}) where B
    i = index(p)
    s.bits[i] && return 0
    s.bits[i] = true

    1 - 2(count(s.bits[1:i-1]) % 2)
end
create(s, p) = (s2 = deepcopy(s); (create!(s2, p), s2))

annihil!(s::Slater{B}, p::Wrapped) where B = annihil!(s, inner(p))
function annihil!(s::Slater{B}, p::MaybeIndex{B}) where B
    i = index(p)
    ~s.bits[i] && return 0
    s.bits[i] = false

    1 - 2(count(s.bits[1:i-1]) % 2)
end
annihil(s, p) = (s2 = deepcopy(s); (annihil!(s2, p), s2))

function Base.show(io::IO, x::Slater)
    SB = innertype(x)
    SBI = indextype(SB)
    str = prod((SBI[I] in x ? string(SB[I]) : "_") * (I == lastindex(SB) ? "" : ", ")
               for I in eachindex(SB))
    print(io, "Slater($(dim(SB)))[$str]")
end
