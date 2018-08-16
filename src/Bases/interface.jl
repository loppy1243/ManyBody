### Interface that an AbstractBasis must implement
# Unexported, optional
#noexport Generation
# Defined in this file
export basis, basistype, dim,  index, indexbasis, inner, innertype

abstract type Generation end
struct Provided <: Generation end
struct Generated <: Generation end
struct Computed <: Generation end
struct Unknown <: Generation end
Generation(::Type{<:AbstractBasis}) = Computed()

basis(B::Type{<:AbstractBasis}) = _basis(B, Generation(B))
_basis(B::Type{<:AbstractBasis}, ::Computed) = B[1:dim(B)]
## Don't know if I want to do this weird stuff with indexbasis...
let _BASIS_CACHE = Dict{Type, Vector{<:AbstractBasis}}()
    global _basis#, _indexbasis
    _basis(B::Type{<:AbstractBasis}, ::Generated) = if haskey(_BASIS_CACHE, B)
        _BASIS_CACHE[B]
    else
        _BASIS_CACHE[B] = B[1:dim(B)]
#        @eval indexbasis(::Type{$B}, i::Int) = _indexbasis($B, i, Generation($B))
    end

#    _indexbasis(B::Type{<:AbstractBasis}, i::Int, ::Generated) = basis(B)[i]
end

index(b::AbstractBasis) = _index(b, Generation(typeof(b)))
_index(b::AbstractBasis, ::Provided) = findfirst(==(b), basis(typeof(b)))

indexbasis(B::Type{<:AbstractBasis}, i::Int) = _indexbasis(B, i, Generation(B))
_indexbasis(B::Type{<:AbstractBasis}, i::Int, ::Provided) = basis(B)[i]

## Do I really want to do this
dim(B::Type{<:AbstractBasis}) = _dim(B, Generation(B))
_dim(B::Type{<:AbstractBasis}, ::Provided) = length(basis(B))

## Really meant as a method on AbstractState's
basistype(x::AbstractBasis) = basistype(typeof(x))
basistype(B::Type{<:AbstractBasis}) = B

## Should be merged with basistype?
innertype(x::AbstractBasis) = innertype(typeof(x))
innertype(B::Type{<:AbstractBasis}) = B

inner(x::AbstractState) = x
