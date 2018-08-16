Base.getindex(::Type{B}, i) where B<:AbstractBasis = _getindex(B, Generation(B), i)
_getindex(::Type{B}, ::Computed, i) where B<:AbstractBasis = indexbasis(B, i)
_getindex(::Type{B}, ::Computed, ixs::Array) where B<:AbstractBasis =
    map(i -> indexbasis(B, i), ixs)
_getindex(::Type{B}, ::Computed, r::AbstractRange{Int}) where B<:AbstractBasis =
    map(i -> indexbasis(B, i), r)

Base.firstindex(::Type{B}) where B<:AbstractBasis = 1
Base.lastindex(::Type{B}) where B<:AbstractBasis = dim(B)

Base.iterate(::Type{B}, i=1) where B<:AbstractBasis = i > dim(B) ? nothing : (B[i], i+1)
Base.IteratorSize(::Type{B}) where B<:AbstractBasis = Base.HasLength()
Base.IteratorEltype(::Type{B}) where B<:AbstractBasis = Base.HasEltype()

Base.length(::Type{B}) where B<:AbstractBasis = dim(B)
Base.eltype(::Type{B}) where B<:AbstractBasis = B
