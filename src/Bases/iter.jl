Base.IteratorSize(::Type{B}) where B<:Type{<:TensorBasis} = Base.HasShape{rank(B)}()
Base.IteratorSize(::Type{B}) where B<:Type{<:AbstractBasis} = Base.HasLength()
Base.axes(B::Type{<:TensorBasis}) = map(d -> 1:d, fulldims(B))
Base.length(B::Type{<:TensorBasis}) = prod(fulldims(B))
Base.length(B::Type{<:AbstractBasis}) = dim(B)

Base.IteratorEltype(::Type{<:Type{<:AbstractBasis}}) = Base.HasEltype()
Base.eltype(B::Type{<:AbstractBasis}) = B

function Base.iterate(B::Type{<:AbstractBasis}, st=(eachindex(B),))
    (x = iterate(st...)) === nothing && return nothing
    (I, inner_st) = x

    (indexer(B)[I], (st[1], inner_st))
end

function Base.iterate(BI::Indexer, st=(eachindex(BI),))
    (x = iterate(st...)) === nothing && return nothing
    (i, inner_st) = x
    (i, (st[1], inner_st))
end
