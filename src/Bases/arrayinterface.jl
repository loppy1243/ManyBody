# If `IndexStyle(Basis) == IndexLinear()`:
#     Must implement
#           `index(::Basis)::Int`
#           `indexbasis(::Type{Basis}, ::Int)::Basis`
# If `IndexStyle(Basis) == IndexCartesian()`:
#       Must implement
#           `index(::Basis)::CartesianIndex`
#           `indexbasis(::Type{Basis}, ::Vararg{Int})::Basis`
# By default, `IndexStyle(::Type{<:TensorBasis}) == IndexCartesian()`.

struct Flat{B<:AbstractBasis} <: AbstractVector{B} end
struct Shaped{N, B<:TensorBasis{N}} <: AbstractArray{B, N}
    Shaped{B}() where B<:TensorBasis = new{rank(B), B}()
end

Flat(B::Type{<:AbstractBasis}) = Flat{B}()
Flat(a::Shaped) = Flat{eltype(a)}()
Shaped(B::Type{<:TensorBasis}) = Shaped{B}()
Shaped(a::Flat{<:TensorBasis}) = Shaped{eltype(a)}()

Base.size(a::Flat) = (dim(eltype(a)),)
Base.size(a::Shaped) = fulldims(eltype(a))

elems(B::Type{<:AbstractBasis}) = Flat(B)
elems(B::Type{<:TensorBasis}) = Shaped(B)

Base.eachindex(B::Type{<:AbstractBasis}) = eachindex(elems(B))
Base.vec(a::Shaped) = Flat(a)
Base.vec(a::Flat) = a
Base.:+(B::Type{<:AbstractBasis}) = elems(B)
Base.:-(B::Type{<:AbstractBasis}) = vec(elems(B))
Base.:+(b::AbstractBasis) = index(b)
Base.:-(b::AbstractBasis) = _minus(b, IndexStyle(b))
_minus(b::AbstractBasis, ::IndexLinear) = index(b)
_minus(b::AbstractBasis, ::IndexCartesian) = LinearIndices(elems(typeof(b)))[index(b)]

indexbasis(B::Type{<:TensorBasis{N}}, I::NTuple{N, Int}) where N = indexbasis(B, I...)
indexbasis(B::Type{<:TensorBasis{N}}, I::CartesianIndex{N}) where N =
    indexbasis(B, Tuple(I)...)

## Could possibly be more efficient to index a TensorBasis linearly, so we account for that.
Base.IndexStyle(::Type{<:AbstractBasis}) = IndexLinear()
Base.IndexStyle(::Type{<:TensorBasis}) = IndexCartesian()

Base.IndexStyle(::Type{<:Flat}) = IndexLinear()
Base.IndexStyle(S::Type{<:Shaped}) = IndexStyle(eltype(S))
Base.IndexStyle(b::AbstractBasis) = IndexStyle(typeof(b))

Base.getindex(a::Flat, i::Int) = indexbasis(eltype(a), i)
Base.getindex(a::Flat{<:TensorBasis}, i::Int) = _getindex(a, i, IndexStyle(eltype(a)))
Base.getindex(a::Shaped, I::Int...) = _getindex(a, I, IndexStyle(a))
_getindex(a::Flat{<:TensorBasis}, i::Int, ::IndexLinear) = indexbasis(eltype(a), i)
_getindex(a::Flat{<:TensorBasis}, i::Int, ::IndexCartesian) =
    indexbasis(eltype(a), CartesianIndices(Shaped(a))[i])
_getindex(a::Shaped, I::Tuple{Int}, ::IndexLinear) = indexbasis(eltype(a), I[1])
_getindex(a::Shaped{N}, I::NTuple{N, Int}, ::IndexCartesian) where N =
    indexbasis(eltype(a), I)

const _Arrayish{B<:AbstractBasis} = Union{Flat{B}, Shaped{<:Any, B}}
Base.in(b::AbstractBasis, a::_Arrayish) = false
Base.in(b::B, a::_Arrayish{B}) where B<:AbstractBasis = true
