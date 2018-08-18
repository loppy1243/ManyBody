struct Index{B<:AbstractBasis} <: AbstractBasis
    index::Int
end
const MaybeIndex{B<:AbstractBasis} = Union{B, Index{B}}
Index(s::AbstractBasis) = Index{typeof(s)}(index(s))

## TODO: Optimize based off of Generated()
Base.:(==)(a::MaybeIndex{B}, b::MaybeIndex{B}) where B<:AbstractBasis = index(a) == index(b)
Base.promote_rule(::Type{Index{B}}, ::Type{B}) where B<:AbstractBasis = B

inner(i::Index) = innertype(i)[i.index]

index(s::Index) = s.index
indexbasis(::Type{S}, s::Int) where S<:Index = S(s)
dim(B::Type{<:Index}) = dim(innertype(B))

Base.convert(::Type{<:Index}, s::AbstractBasis) = Index(s)
#Base.convert(::Type{B}, s::Index{B}) where B<:AbstractBasis = B(s)
Base.convert(::Type{C}, s::Index{B}) where {C<:AbstractBasis, B<:C} = B(s)

B(s::Index{B}) where B<:AbstractBasis = B[s.index]
