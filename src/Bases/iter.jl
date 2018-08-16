Base.getindex(B::Type{<:AbstractBasis}, i) = indexbasis(B, i)
Base.getindex(B::Type{<:AbstractBasis}, ixs::Array) =
   map(i -> indexbasis(B, i), ixs)
Base.getindex(B::Type{<:AbstractBasis}, r::AbstractRange{Int}) =
   map(i -> indexbasis(B, i), r)

Base.firstindex(B::Type{<:AbstractBasis}) = 1
Base.lastindex(B::Type{<:AbstractBasis}) = dim(B)

Base.iterate(B::Type{<:AbstractBasis}, i=1) = i > dim(B) ? nothing : (B[i], i+1)
Base.IteratorSize(B::Type{<:AbstractBasis}) = Base.HasLength()
Base.IteratorEltype(B::Type{<:AbstractBasis}) = Base.HasEltype()

Base.length(B::Type{<:AbstractBasis}) = dim(B)
Base.eltype(B::Type{<:AbstractBasis}) = B
