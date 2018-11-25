Base.IteratorSize(::Type{B}) where B<:Type{<:TensorBasis} = Base.HasShape{rank(B)}()
Base.IteratorSize(::Type{B}) where B<:Type{<:AbstractBasis} = Base.HasLength()
Base.axes(B::Type{<:TensorBasis}) = map(d -> 1:d, fulldims(B))
Base.length(B::Type{<:TensorBasis}) = prod(fulldims(B))
Base.length(B::Type{<:AbstractBasis}) = dim(B)

Base.size(B::Type{<:TensorBasis}) = fulldims(B)
Base.size(B::Type{<:AbstractBasis}) = (dim(B),)

Base.IteratorEltype(::Type{<:Type{<:AbstractBasis}}) = Base.HasEltype()
Base.eltype(B::Type{<:AbstractBasis}) = B
#Base.eltype(B::Type{<:Sub}) = supbasis(B)
Base.eltype(B::Type{<:Indexer{SB}}) where SB<:Sub = supbasis(SB)
Base.eltype(B::Type{<:LinearIndexer{SB}}) where SB<:Sub = supbasis(SB)

function Base.iterate(B::Type{<:AbstractBasis}, st=(eachindex(B),))
    (x = iterate(st...)) === nothing && return nothing
    (I, inner_st) = x

    (indexer(B)[I], (st[1], inner_st))
end
