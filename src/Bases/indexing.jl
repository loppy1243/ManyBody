Base.getindex(B::Type{<:AbstractBasis}, ixs...) = indexbasis(B, ixs...)
Base.to_index(b::AbstractBasis) = index(b)

Base.eachindex(B::Type{<:TensorBasis}) = CartesianIndices(axes(B))
Base.firstindex(B::Type{<:TensorBasis}) = CartesianIndex(ones(Int, rank(B))...)
Base.lastindex(B::Type{<:TensorBasis}) = CartesianIndex(fulldims(B))

Base.CartesianIndices(B::Type{<:AbstractBasis}) = CartesianIndices((dim(B),))
Base.LinearIndices(B::Type{<:AbstractBasis}) = LinearIndices((dim(B),))
Base.CartesianIndices(B::Type{<:TensorBasis}) = CartesianIndices(fulldims(B))
Base.LinearIndices(B::Type{<:TensorBasis}) = LinearIndices(fulldims(B))

linearindex(B::Type{<:AbstractIndex}) = LinearIndices(b)[b]
