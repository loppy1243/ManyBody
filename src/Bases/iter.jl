Base.getindex(B::Type{<:AbstractBasis}, ixs...) = indexbasis(B, ixs...)
Base.getindex(B::Type{<:AbstractBasis}, ixs::AbstractArray) =
   map(i -> indexbasis(B, i), ixs)

function Base.eachindex(B::Type{<:ConcreteBasis})
    _eachindex(::IndexTypes.Linear) = 1:dim(B)
    _eachindex(::IndexTypes.Cartesian) = Base.CartesianIndices(innerdims(B))

    _eachindex(IndexType(B))
end

Base.firstindex(B::Type{<:AbstractBasis}) = first(eachindex(B))
Base.lastindex(B::Type{<:AbstractBasis}) = last(eachindex(B))

function Base.iterate(B::Type{<:AbstractBasis})
    ix_iter = eachindex(B)
    (ix_iter_ret = iterate(ix_iter)) === nothing && return nothing

    (B[ix_iter_ret[1]], (ix_iter, ix_iter_ret[2]))
end
function Base.iterate(B::Type{<:AbstractBasis}, st)
    ix_iter, ix_iter_st = st
    (ix_iter_ret = iterate(ix_iter, ix_iter_st)) === nothing && return nothing

    (B[ix_iter_ret[1]], (ix_iter, ix_iter_ret[2]))
end

Base.IteratorSize(::Type{<:Type{<:AbstractBasis}}) = Base.HasLength()
Base.IteratorEltype(::Type{<:Type{<:AbstractBasis}}) = Base.HasEltype()

Base.length(B::Type{<:AbstractBasis}) = dim(B)
Base.eltype(B::Type{<:AbstractBasis}) = B
