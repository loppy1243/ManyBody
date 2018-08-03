export indexes

struct Index{B<:AbstractBasis} <: AbstractBasis
    index::Int
end
Index(s::AbstractBasis) = Index{typeof(s)}(index(s))

indexes(::Type{B}) where B<:AbstractBasis = map(Index{B}, indices(B))

index(s::Index) = s.index
indexbasis(::Type{S}, s::Int) where S<:Index = S(s)
basis(::Type{Index{B}}) where B = indexes(B)
dim(::Type{Index{B}}) where B = dim(B)

Base.convert(::Type{Index{B}}, s::B) where B = Index{B}(index(s))
Base.convert(::Type{Index}, s::AbstractBasis) = Index{typeof(s)}(index(s))
Base.convert(::Type{B}, s::Index{B}) where B<:AbstractBasis = B(s)

B(s::Index{B}) where B<:AbstractBasis = B[s.index]
