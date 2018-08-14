### Interface that an AbstractBasis must implement
# Unexported, optional
#noexport Generation
# Defined in this file, optional
export basis
# Not defined here
export index, indexbasis, dim

abstract type Generation end
struct Provided <: Generation end
struct Generated <: Generation end
struct Computed <: Generation end
Generation(::Type{<:AbstractBasis}) = Computed()

basis(::Type{B}) where B<:AbstractBasis = _basis(B, Generation(B))
_basis(::Type{B}, ::Computed) where B<:AbstractBasis = B[1:dim(B)]
@generated _basis(::Type{B}, ::Generated) where B<:AbstractBasis =
    map(i -> indexbasis(B, i), 1:dim(B))
