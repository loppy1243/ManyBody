struct Index{B<:AbstractBasis} <: AbstractBasis
    index::Int
end
const MaybeIndex{B<:AbstractBasis} = Union{B, Index{B}}
Index(s::AbstractBasis) = Index{typeof(s)}(index(s))

Base.:(==)(a::MaybeIndex{B}, b::MaybeIndex{B}) where B<:AbstractBasis = index(a) == index(b)
Base.promote_rule(::Type{Index{B}}, ::Type{B}) where B<:AbstractBasis = B

inner(i::Index) = innertype(i)[i.index]

index(s::Index) = s.index
indexbasis(S::Type{<:Index}, s::Int) = S(s)
dim(B::Type{<:Index}) = dim(innertype(B))

Base.convert(::Type{Index}, s::AbstractBasis) = Index(s)
Base.convert(::Type{Index{B}}, s::B) where B<:AbstractBasis = Index{B}(s)
Base.convert(I::Type{<:Index}, s::AbstractBasis) = Index(convert(innertype(I), s))
Base.convert(::Type{B}, s::Index{B}) where B<:AbstractBasis = B[index(s)]
