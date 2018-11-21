using Base.Cartesian: @ncall

struct Indexer{B<:AbstractBasis, N} <: AbstractArray{B, N} end
Indexer{B}() where B<:TensorBasis = Indexer{B, rank(B)}()
Indexer{B}() where B<:AbstractBasis = Indexer{B, 1}()

struct LinearIndexer{B<:AbstractBasis} <: AbstractVector{B} end

indexer(B::Type{<:AbstractBasis}) = Indexer{B}()
linearindexer(B::Type{<:AbstractBasis}) = LinearIndexer{B}()
_basis(::Type{<:Indexer{B}}) where B<:AbstractBasis = B
_basis(BI::Indexer) = _basis(typeof(BI))
_basis(::Type{<:LinearIndexer{B}}) where B<:AbstractBasis = B
_basis(BI::LinearIndexer) = _basis(typeof(BI))

## Fed up with this shit...
Base.show(io::IO, x::Indexer) = show(io, collect(x))
Base.show(io::IO, mime::MIME"text/plain", x::Indexer) = show(io, mime, collect(x))

indexbasis(BI::Indexer, ixs::Int...) = indexbasis(_basis(BI), ixs...)
indexbasis(BI::LinearIndexer, ix::Int) = indexbasis(_basis(BI), ix)
indexbasis(BI::LinearIndexer{<:TensorBasis}, ix::Int) =
    indexbasis(_basis(BI), Tuple(CartesianIndices(_basis(BI))[ix])...)

Base.size(BI::Indexer) = (dim(_basis(BI)),)
Base.size(BI::Indexer{<:TensorBasis}) = fulldims(_basis(BI))
Base.size(BI::LinearIndexer) = (dim(_basis(BI)),)
IndexStyle(::Type{<:Indexer}) = Base.IndexLinear()
IndexStyle(::Type{<:Indexer{<:TensorBasis}}) = Base.IndexCartesian()
IndexStyle(::Type{<:LinearIndexer}) = Base.IndexLinear()

Base.getindex(BI::Indexer, ixs...) = indexbasis(_basis(BI), Base.to_indices(BI, ixs)...) 
Base.getindex(BI::LinearIndexer, ix) =
    indexbasis(_basis(BI),
               Tuple(CartesianIndices(_basis(BI))[Base.to_index(BI, ix)...])...)

Base.to_index(::Any, b::AbstractBasis) = index(b)
Base.to_indices(A, x::Tuple{<:TensorBasis, Vararg{Any}}) =
    (Tuple(index(x[1]))..., Base.to_indices(A, Base.tail(x))...)

Base.eachindex(BI::Indexer) = eachindex(_basis(BI))
Base.eachindex(BI::LinearIndexer) = 1:dim(_basis(BI))
Base.eachindex(B::Type{<:AbstractBasis}) = LinearIndices(B)
Base.eachindex(B::Type{<:TensorBasis}) = CartesianIndices(B)

Base.firstindex(BI::Indexer) = firstindex(_basis(BI))
Base.firstindex(::LinearIndexer) = 1
Base.firstindex(B::Type{<:AbstractBasis}) = 1
Base.firstindex(B::Type{<:TensorBasis}) = CartesianIndex(ones(Int, rank(B))...)

Base.lastindex(BI::Indexer) = lastindex(_basis(BI))
Base.lastindex(BI::LinearIndexer) = dim(_basis(BI))
Base.lastindex(B::Type{<:AbstractBasis}) = dim(B)
Base.lastindex(B::Type{<:TensorBasis}) = CartesianIndex(fulldims(B))

Base.CartesianIndices(BI::Indexer) = CartesianIndices(_basis(BI))
Base.CartesianIndices(BI::LinearIndexer) = CartesianIndices((dim(_basis(BI)),))
Base.CartesianIndices(B::Type{<:AbstractBasis}) = CartesianIndices((dim(B),))
Base.CartesianIndices(B::Type{<:TensorBasis}) = CartesianIndices(axes(B))

Base.LinearIndices(BI::Indexer) = LinearIndices(_basis(BI))
Base.LinearIndices(BI::LinearIndexer) = LinearIndices(_basis(BI))
Base.LinearIndices(B::Type{<:AbstractBasis}) = LinearIndices((dim(B),))
Base.LinearIndices(B::Type{<:TensorBasis}) = LinearIndices(axes(B))

"""
    linearindex(b)

Compute the linear index for the basis element `b`.

Usually you want to use `LinearIndices(typeof(b))` instead so that the linear index of
multiple elements may be computed with less overhead.
"""
linearindex(b::AbstractBasis) = LinearIndices(typeof(b))[b]
