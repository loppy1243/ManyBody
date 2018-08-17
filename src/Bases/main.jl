@reexport module Bases
export RefStates, AbstractState, AbstractBasis, RefState, Bra, ZeroState

using ..AbstractState
using Base.Cartesian: @ncall
using LinearAlgebra: Adjoint
using Combinatorics: combinations
using Reexport: @reexport

abstract type AbstractBasis <: AbstractState end

struct Bra{S}; state::S end
struct ZeroState <: AbstractState end

include("interface.jl")
include("iter.jl")
include("indexbasis.jl")
include("subbasis.jl")

## FIXME: Move Neg so that this works
const Rep{B<:AbstractBasis} = Union{B, Sub{B}, Index{B}, Neg{B}, Product{1, Tuple{B}}}

include("pairing.jl")
include("refstates.jl")
include("slater.jl")
include("product.jl")
include("baseops/main.jl")

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

end # module States
