export indexes

struct Index{SP<:SPBasis} <: SPBasis
    index::Int
end
Index(s::SPBasis) = Index{typeof(s)}(index(s))

indexes(::Type{B}) where B<:Basis = map(Index{B}, indices(B))

index(s::Index) = s.index
indexbasis(::Type{S}, s::Int) where S<:Index = S(s)
basis(::Type{Index{B}}) where B = indexes(B)
dim(::Type{Index{B}}) where B = dim(B)

Base.convert(::Type{Index{B}}, s::B) where B = Index{B}(index(s))
Base.convert(::Type{Index}, s::Basis) = Index{typeof(s)}(index(s))
Base.convert(::Type{B}, s::Index{B}) where B<:Basis = B(s)

SP(s::Index{SP}) where SP<:SPBasis = SP[s.index]
