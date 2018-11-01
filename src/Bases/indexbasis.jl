export indextype, rank
using ..ManyBody: @disallow

module IndexTypes
    export IndexType

    abstract type IndexType end
    struct Cartesian{N} <: IndexType end
    struct Linear <: IndexType end

    rank(::Linear) = 1
    rank(::Cartesian{N}) where N = N
end
using .IndexTypes

rank(b) = rank(typeof(b))
@disallow rank(::Type)
rank(B::Type{<:ConcreteBasis}) = IndexTypes.rank(IndexType(B))

indextype(b) = indextype(typeof(b))
@disallow indextype(::Type)
indextype(B::Type{<:ConcreteBasis}) = indextype(B, IndexType(B))
indextype(B::Type{<:ConcreteBasis}, ::IndexTypes.Linear) = LinearIndex{B}
indextype(B::Type{<:ConcreteBasis}, ::IndexTypes.Cartesian) = CartesianIndex{B, rank(B)}

struct LinearIndex{B<:ConcreteBasis} <: AbstractIndex{B}
    index::Int

    LinearIndex{B}(index::Int, ::IndexTypes.Linear) where B<:ConcreteBasis = new(index)
end
LinearIndex{B}(index::Int) where B<:ConcreteBasis = LinearIndex{B}(index, IndexType(B))
LinearIndex{B}(index::Tuple{Int}) where B<:ConcreteBasis =
    LinearIndex{B}(index[1], IndexType(B))

struct CartesianIndex{B<:ConcreteBasis, N} <: AbstractIndex{B}
    index::Base.CartesianIndex{N}

    CartesianIndex{B, N}(index, ::IndexTypes.Cartesian{N}) where {B<:ConcreteBasis, N} =
        new(index)
end
CartesianIndex{B, N}(index) where {B<:ConcreteBasis, N} =
    CartesianIndex{B, N}(index, IndexType(B))
CartesianIndex{B}(index) where B<:ConcreteBasis = CartesianIndex{B, length(index)}(index)
CartesianIndex{B}(indices::Int...) where B<:ConcreteBasis =
    CartesianIndex{B, length(indices)}(indices)

dim(I::Type{<:AbstractIndex}) = dim(innertype(I))
innerdims(I::Type{AbstractIndex}) = innerdims(innertype(I))

basis(I::Type{<:LinearIndex}) = map(I, eachindex(basistype(I)))

index(I::LinearIndex) = I.index
index(I::CartesianIndex) = I.index

indexbasis(::Type{B}, I::AbstractIndex{B}) where B<:ConcreteBasis = indexbasis(B, index(I))
indexbasis(I::Type{<:AbstractIndex}, ixs...) = I(ixs...)
indexbasis(B::Type{<:AbstractBasis}, I::AbstractIndex) = convert(B, I)

innertype(::Type{<:AbstractIndex{B}}) where B<:ConcreteBasis = B
inner(I::AbstractIndex) = innertype(I)[I]

Base.:(==)(a::I, b::I) where I<:AbstractIndex = index(a) == index(b)

Base.convert(B::Type{<:ConcreteBasis}, I::AbstractIndex) = convert(B, inner(I))
Base.convert(I::Type{<:AbstractIndex}, b::AbstractIndex) = convert(I, convert(innertype(I), inner(b)))

Base.convert(I::Type{<:AbstractIndex}, b::ConcreteBasis) = I(index(convert(innertype(I), b)))
Base.convert(::Type{AbstractIndex}, b::ConcreteBasis) = indextype(B)(index(b))
Base.convert(::Type{LinearIndex{B}}, b::B) where B<:ConcreteBasis = LinearIndex{B}(index(b))
Base.convert(::Type{<:CartesianIndex{B}}, b::B) where B<:ConcreteBasis =
    CartesianIndex{B}(index(b))
#Base.convert(::Type{LinearIndex{B}}, b::LinearIndex{B}) where B<:ConcreteBasis = b 

## These methods seems to be missing from Base.
## Causes StackOverflowError?
#Base.convert(::Type{<:Base.CartesianIndex{1}}, x::Integer) = Base.CartesianIndex(x)
#Base.convert(::Type{<:Base.CartesianIndex{N}}, x::NTuple{N, Integer}) where N =
#    Base.CartesianIndex(x)

#Base.convert(I::Type{<:LinearIndex}, x) = I(convert(Int, x))
#Base.convert(I::Type{<:CartesianIndex}, x) = I(convert(Base.CartesianIndex, x))
