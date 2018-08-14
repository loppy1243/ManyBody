export create, create!, annihil, annihil!

struct Slater{B<:AbstractBasis} <: AbstractBasis
    bits::BitVector

    Slater{B}(parts::BitVector) where B = new(parts)
    function Slater{B}(parts) where B
        bits = falses(dim(B))

        for p in parts
            p = convert(B, p)

            bits[index(p)] = true
        end

        new(bits)
    end
end
Slater(x::Tuple) = MethodError(Slater, (x,)) |> throw
Slater(ps...) = Slater(ps)
Slater(ps::NTuple{N, B}) where {N, B<:AbstractBasis} = Slater{B}(ps)

Base.:(==)(s1::MB, s2::MB) where MB<:Slater = s1.bits == s2.bits
Base.in(p::B, s::Slater{B}) where B = s.bits[index(p)]
Base.in(p::B, s::Sub{Slater{B}}) where B = p in s.state

index(s::Slater) = sum(s.bits[i]*2^(i-1) for i in eachindex(s.bits))
function indexbasis(::Type{Slater{B}}, ix::Int) where B 
    d = dim(B)
    states = falses(d)
    for i in 1:d
        states[i] = (ix & (1 << (i-1))) != 0
    end

    Slater{B}(states)
end

dim(::Type{Slater{B}}) where B = 2^dim(B) - 1
n_occ(s::Slater) = count(s.bits)
n_unoc(s::Slater) = count(.~s.bits)

function create!(s::Slater{B}, p::B) where B
    i = index(p)
    s.bits[i] = true

    1 - 2(count(s.bits[1:i-1]) % 2)
end
function annihil!(s::Slater{B}, p::B) where B
    i = index(p)
    s.bits[i] = false

    1 - 2(count(s.bits[1:i-1]) % 2)
end
create(s, p) = (s2 = deepcopy(s); (create!(s2, p), s2))
annihil(s, p) = (s2 = deepcopy(s); (annihil!(s2, p), s2))

function Base.show(io::IO, x::MaybeSub{Slater{B}}) where B
    x = convert(Slater{B}, x)
    print(io, "Slater($(dim(B)))[",
              [string(B[i])*" " for (i, p) in enumerate(x.bits) if p]...,
              "]")
end