@reexport module Bases
#export AbstractBasis, ConcreteBasis, LeafBasis, RefStates, RefState

import ..ManyBody
using ..ManyBody: AbstractState
using Base.Cartesian: @ncall
using LinearAlgebra: Adjoint
using Combinatorics: combinations
using Reexport: @reexport

#abstract type AbstractBasis <: AbstractState end
#abstract type ConcreteBasis <: AbstractBasis end
#abstract type LeafBasis <: AbstractBasis end
#abstract type AbstractIndex{B<:ConcreteBasis} <: LeafBasis end
#const MaybeIndex{B<:ConcreteBasis} = Union{B, AbstractIndex{B}}

abstract type AbstractBasis
abstract type TensorBasis{Rank} <: AbstractBasis end
const Basis = TensorBasis{1}
#abstract type FockBasis{SPBasis<:TensorBasis} <: AbstractBasis end

rank(::Type{<:TensorBasis{N}}) where N = N
rank(b::TensorBasis) = rank(typeof(b))
dim(B::Type{<:TensorBasis}) = prod(fulldims(B))
# Good way to do this?
#fulldims(B::Type{<:Basis}) = (dim(B),)

Base.:(==)(x::B, y::B) where B<:AbstractBasis = index(x) == index(y)

Base.Array{T}(b::TensorBasis) where T = Array{T, rank(b)}(b)
function Base.Array{T, N}(b::TensorBasis{N}) where {T, N}
    ret = zeros(T, fulldims(typeof(b)))
    ret[b] = oneunit(T)

    ret
end
function Base.Vector{T}(b::AbstractBasis) where T
    ret = zeros(T, dim(typeof(b)))
    ret[linearindex[b]] = oneunit(T)

    ret
end

include("iter.jl")
include("pairing.jl")
include("subbasis.jl")
include("product.jl")
include("slater.jl")

### Update Line
##############################################################################################

include("refstates.jl")
end # module States
