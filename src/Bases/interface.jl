### Interface that an AbstractBasis must implement
# Unexported, optional
#noexport Generation
# Defined in this file
export basis, basistype, dim
# Not defined here
export index, indexbasis

abstract type Generation end
struct Provided <: Generation end
struct Generated <: Generation end
struct Computed <: Generation end
struct Unknown <: Generation end
Generation(::Type{<:AbstractBasis}) = Computed()

basis(B::Type{<:AbstractBasis}) = _basis(B, Generation(B))
_basis(B::Type{<:AbstractBasis}, ::Computed) = B[1:dim(B)]
let _BASIS_CACHE = Dict{Type, Vector{<:AbstractBasis}}()
    global _basis(::Type{B}, ::Generated) where B<:AbstractBasis = if haskey(_BASIS_CACHE, B)
        _BASIS_CACHE[B]
    else
        _BASIS_CACHE[B] = map(i -> indexbasis(B.parameters[1], i), 1:dim(B))
    end
end

dim(B::Type{<:AbstractBasis}) = _dim(B, Generation(B))
_dim(::Type{B}, ::Provided) where B<:AbstractBasis = length(basis(B))

basistype(x::AbstractBasis) = basistype(typeof(x))
basistype(::Type{B}) where B<:AbstractBasis = B
