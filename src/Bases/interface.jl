### Interface that an AbstractBasis must implement
export basis, innertype, innerypes, inner, innerdims, constrain
#export dim, index, indexbasis

## FIXME: Move Neg so that this works
const Wrapped{B<:ConcreteBasis} = Union{Sub{B}, AbstractIndex{B}, Neg{B}, Product{1, Tuple{B}}}
const Rep{B<:ConcreteBasis} = Union{B, Wrapped{B}}

#Base.:(==)(a::Rep{B}, b::Rep{B}) where B<:ConcreteBasis = inner(a) == inner(b)

basis(B::Type{<:AbstractBasis}) = B[eachindex(B)]

innertype(x::AbstractBasis) = innertype(typeof(x))
innertype(B::Type{<:AbstractBasis}) = B
innertypes(B::Type{<:AbstractBasis}) = (innertype(B),)

innerdims(B::Type{<:AbstractBasis}) = (dim(B),)

inner(x::AbstractBasis) = x

Base.convert(B1::Type{<:ConcreteBasis}, b2::Wrapped) = Base.convert(B1, inner(b2))

constrain(::Type{B}, b::B) where B<:ConcreteBasis = b
constrain(SB::Type{<:Sub}, b::ConcreteBasis) = SB(constrain(innertype(SB), b))

constrain(B::Type{<:ConcreteBasis}, b::AbstractIndex) = indextype(B)[b]
constrain(B, b::Product{1, Tuple{B2}}) where B2<:ConcreteBasis = constrain(B, inner(b))
constrain(B, b::Neg) = Neg(constrain(B, inner(b)))
