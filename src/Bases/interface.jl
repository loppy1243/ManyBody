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
@generated _basis(::Type{B}, ::Generated) where B<:AbstractBasis =
    map(i -> indexbasis(B.parameters[1], i), 1:dim(B))

dim(B::Type{<:AbstractBasis}) = _dim(B, Generation(B))
_dim(::Type{B}, ::Provided) where B<:AbstractBasis = length(basis(B))

basistype(x::AbstractBasis) = basistype(typeof(x))
basistype(::Type{B}) where B<:AbstractBasis = B
