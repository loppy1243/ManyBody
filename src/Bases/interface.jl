export basis, innertype, innertypes, inner, innerdims, widen, constrain
## Not defined here
export dim, index, indexbasis

const Wrapped{B<:ConcreteBasis} = Union{Sub{B}, AbstractIndex{B}, Neg{B}, Product{1, Tuple{B}}}
const Rep{B<:ConcreteBasis} = Union{B, Wrapped{B}}

#Base.:(==)(a::Rep{B}, b::Rep{B}) where B<:ConcreteBasis = inner(a) == inner(b)

basis(B::Type{<:AbstractBasis}) = B[eachindex(B)]

innertype(x::AbstractBasis) = innertype(typeof(x))
innertype(B::Type{<:AbstractBasis}) = B
innertypes(B::Type{<:AbstractBasis}) = (innertype(B),)

innerdims(B::Type{<:AbstractBasis}) = (dim(B),)

inner(x::AbstractBasis) = x

#Base.convert(B1::Type{<:ConcreteBasis}, b2::Wrapped) = Base.convert(B1, inner(b2))

widen(::Type{B}, b::B) where B<:ConcreteBasis = b
widen(B::Type{<:ConcreteBasis}, sb::Sub) = Bases.widen(B, inner(sb))
widen(B::Type{<:Product{N}}, sb::Product{N}) where N =
    Product(map(Bases.widen, innertypes(B), sb.states))

widen(B::Type{<:ConcreteBasis}, b::AbstractIndex) =
    convert(indextype(B), Bases.widen(B, inner(b)))
widen(B::Type{<:Product{1}}, b::ConcreteBasis) = Product(Bases.widen(innertype(B), b))
widen(B::Type{<:ConcreteBasis}, b::Neg) = Neg(Bases.widen(B, inner(b)))

constrain(::Type{B}, b::B) where B<:ConcreteBasis = b
constrain(SB::Type{<:Sub{B}}, b::B) where B<:ConcreteBasis = SB(b)
constrain(SB::Type{<:Sub}, b::ConcreteBasis) = SB(constrain(innertype(SB), b))
constrain(SB::Type{<:Product}, b::Product) = Product(map(constrain, innertypes(SB), b.states))

constrain(B::Type{<:ConcreteBasis}, b::AbstractIndex) =
    convert(indextype(B), (constrain(B, inner(b))))
constrain(B::Type{<:ConcreteBasis}, b::Product{1}) where B2<:ConcreteBasis =
    constrain(B, inner(b))
constrain(B::Type{<:ConcreteBasis}, b::Neg) = Neg(constrain(B, inner(b)))
