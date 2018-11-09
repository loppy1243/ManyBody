Base.getindex(B::Type{<:TensorBasis}, ixs...) = indexbasis(B, ixs...)
Base.to_index(b::TensorBasis) = index(b)

Base.eachindex(B::Type{<:TensorBasis}) = CartesianIndices(axes(B))
Base.firstindex(B::Type{<:TensorBasis}) = CartesianIndex(ones(rank(B))...)
Base.lastindex(B::Type{<:TensorBasis}) = CartesianIndex(fulldims(B))

Base.IteratorSize(::Type{B}) where B<:Type{<:TensorBasis} = Base.HasShape{rank(B)}()
Base.axes(B::Type{<:TensorBasis}) = map(d -> 1:d, fulldims(B))
Base.length(B::Type{<:TensorBasis}) = prod(fulldims(B))

Base.IteratorEltype(::Type{<:Type{<:TensorBasis}}) = Base.HasEltype()
Base.eltype(B::Type{<:TensorBasis}) = B

function Base.iterate(B::Type{<:TensorBasis}, st=(eachindex(B),))
    (x = iterate(st...)) === nothing && return nothing
    (I, inner_st) = x

    (B[I], (itr, inner_st))
end
