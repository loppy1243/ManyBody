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
rank(b::AbstractBasis) = rank(typeof(b))
dim(B::Type{<:TensorBasis}) = prod(fulldims(B))
# Good way to do this?
#fulldims(B::Type{<:Basis}) = (dim(B),)

Base.:(==)(x::B, y::B) where B<:AbstractBasis = index(x) == index(y)

include("iter.jl")
include("pairing.jl")
include("subbasis.jl")
include("product.jl")

### Update Line
##############################################################################################

include("slater.jl")
include("refstates.jl")
end # module States
