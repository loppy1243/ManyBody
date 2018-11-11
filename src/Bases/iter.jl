Base.IteratorSize(::Type{B}) where B<:Type{<:TensorBasis} = Base.HasShape{rank(B)}()
Base.axes(B::Type{<:TensorBasis}) = map(d -> 1:d, fulldims(B))
Base.length(B::Type{<:TensorBasis}) = prod(fulldims(B))

Base.IteratorEltype(::Type{<:Type{<:AbstractBasis}}) = Base.HasEltype()
Base.eltype(B::Type{<:AbstractBasis}) = B

function Base.iterate(B::Type{<:AbstractBasis}, st=(eachindex(B),))
    (x = iterate(st...)) === nothing && return nothing
    (I, inner_st) = x

    (B[I], (itr, inner_st))
end
