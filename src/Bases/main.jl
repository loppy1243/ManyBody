@reexport module Bases
export AbstractBasis, ConcreteBasis, LeafBasis, RefStates, RefState

import ..ManyBody
using ..ManyBody: AbstractState
using Base.Cartesian: @ncall
using LinearAlgebra: Adjoint
using Combinatorics: combinations
using Reexport: @reexport

abstract type AbstractBasis <: AbstractState end
abstract type ConcreteBasis <: AbstractBasis end
abstract type LeafBasis <: AbstractBasis end
abstract type AbstractIndex{B<:ConcreteBasis} <: LeafBasis end
const MaybeIndex{B<:ConcreteBasis} = Union{B, AbstractIndex{B}}

ManyBody.basistype(x::AbstractBasis) = basistype(typeof(x))
ManyBody.basistype(B::Type{<:AbstractBasis}) = B

include("indexbasis.jl")
include("subbasis.jl")
include("product.jl")
include("interface.jl")
include("iter.jl")
include("pairing.jl")
include("refstates.jl")
include("slater.jl")

@defSub Paired{F, L} <: Slater{Pairing{L}} begin
    SP = Pairing{L}
    ref = Slater{SP}(SP(l, s) for l = 1:F, s in SPINS)

    ret = [ref]
    for ph = 1:L-F
        for lhs in combinations(1:F, ph), lps in combinations(F+1:L, ph)
            state = deepcopy(ref)
            for (lh, lp) in zip(lhs, lps), s in SPINS
                annihil!(state, SP(lh, s))
                create!(state, SP(lp, s))
            end

            push!(ret, state)
        end
    end

    ret
end
IndexType(::Type{<:Paired}) = IndexTypes.Linear()

end # module States
