#@generated function Base.getindex(B::Type{<:AbstractBasis}, ixs...)
#    full_ixs = []
#    for (i, T) in enumerate(ixs)
#        if T <: CartesianIndex
#            append!(full_ixs, [:(Tuple(ixs[$i])[$j]) for j = 1:length(T)])
#        else
#            push!(full_ixs, :(ixs[$i]))
#        end
#    end
#
#    :(indexbasis(B, $(full_ixs...)))
#end
function Base.getindex(B::Type{<:AbstractArray}, ixs...)
    N = sum(map(length, ixs))
    full_ixs = Vector{Int}(undef, N)
    for (i, jxs) in enumerate(ixs)
        for (j, jx) in enumerate(Bases.to_indices(B, jxs))
            full_ixs[i+j-1] = jx
        end
    end

    indexbasis(B, full_ixs...)
end
to_indices(::Type{<:AbstractBasis}, x) = (x,)
to_indices(::Type{<:TensorBasis}, x::CartesianIndex) = Tuple(x)

Base.to_indices(::Any, b::AbstractBasis) = (index(b),)
Base.to_indices(::Any, b::TensorBasis) = Tuple(index(b))

Base.eachindex(B::Type{<:AbstractBasis}) = LinearIndices(B)
Base.firstindex(B::Type{<:AbstractBasis}) = 1
Base.lastindex(B::Type{<:AbstractBasis}) = dim(B)

Base.eachindex(B::Type{<:TensorBasis}) = CartesianIndices(B)
Base.firstindex(B::Type{<:TensorBasis}) = CartesianIndex(ones(Int, rank(B))...)
Base.lastindex(B::Type{<:TensorBasis}) = CartesianIndex(fulldims(B))

Base.CartesianIndices(B::Type{<:AbstractBasis}) = CartesianIndices((dim(B),))
Base.LinearIndices(B::Type{<:AbstractBasis}) = LinearIndices((dim(B),))
Base.CartesianIndices(B::Type{<:TensorBasis}) = CartesianIndices(axes(B))
Base.LinearIndices(B::Type{<:TensorBasis}) = LinearIndices(axes(B))

"""
    linearindex(b)

Compute the linear index for the basis element `b`.

Usually you want to use `LinearIndices(typeof(b))` instead so that the linear index of
multiple elements may be computed with less overhead.
"""
linearindex(B::AbstractBasis) = LinearIndices(typeof(b))[b]
