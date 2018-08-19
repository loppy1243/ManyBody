### Interface that an AbstractBasis must implement
export basis, basistype, supbasistype, dim, index, indexbasis, inner, innertype
using ..@disallow

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
basistype(x) = basistype(typeof(x))
@disallow basistype(::Type)
basistype(x::AbstractBasis) = basistype(typeof(x))
basistype(B::Type{<:AbstractBasis}) = B

supbasistype(B::Type{<:AbstractBasis}) = B
supbasistype(B::Type{<:Wrapped}) = supbasistype(innertype(B))
supbasistype(T) = supbasistype(basistype(T))

innertype(x::AbstractBasis) = innertype(typeof(x))
innertype(B::Type{<:AbstractBasis}) = B
innertype(::Type{<:Wrapped{B}}) where B<:AbstractBasis = B

inner(x::AbstractBasis) = x

indextype(B::Type{<:AbstractBasis}) = Index{B}
indextype(I::Type{<:Index}) = I




Base.convert(::Type{B}, b2::Rep{B}) where B<:AbstractBasis = inner(b2)
Base.convert(B1::Type{<:AbstractBasis}, b2::Wrapped) = Base.convert(B1, inner(b2))
