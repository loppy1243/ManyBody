@reexport module Bases
export RefStates, AbstractState, AbstractBasis, RefState, Bra, overlap
# Basis interface
export basis, index, indexbasis, dim #, Base.indices

using ..AbstractState
using Reexport: @reexport

abstract type AbstractBasis <: AbstractState end

Base.:(==)(x::AbstractBasis, y::AbstractBasis) = false
Base.:(==)(x::B, y::B) where B<:AbstractBasis = index(x) == index(y)

Base.getindex(::Type{B}, ixs...) where B<:AbstractBasis = indexbasis(B, ixs...)
Base.getindex(::Type{B}, ixs::Array) where B<:AbstractBasis = map(i -> indexbasis(B, i), ixs)
Base.getindex(::Type{B}, r::Range{Int}) where B<:AbstractBasis = map(i -> indexbasis(B, i), r)

Base.start(::Type{B}) where B<:AbstractBasis = 1
Base.next(::Type{B}, st) where B<:AbstractBasis = (B[st], st+1)
Base.done(::Type{B}, st) where B<:AbstractBasis = st > dim(B)

Base.endof(::Type{B}) where B<:AbstractBasis = dim(B)
Base.indices(::Type{B}) where B<:AbstractBasis = 1:dim(B)
# NOTE: May not be a good idea for large bases
@generated basis(::Type{B}) where B<:AbstractBasis = B[indices(B)]

include("indexbasis.jl")
include("subbasis.jl")
include("pairing.jl")
include("refstates.jl")
include("mbbasis.jl")

@defSubBasis Paired{F, L} <: Slater{Pairing{L}} begin
    SP = Pairing{L}
    ref = Slater{SP}(SP(l, s) for l = 1:F, s in SPINS)

    ret = [ref]
    for ph = 1:L-F
        for lhs in combinations(1:F, ph), lps in combinations(F+1:L, ph)
            state = deepcopy(ref)
            for (lh, lp) in zip(lhs, lps)
                annihil!(state, SP(lh, SPINUP))
                annihil!(state, SP(lh, SPINDOWN))
                create!(state, SP(lp, SPINUP))
                create!(state, SP(lp, SPINUP))
            end

            push!(ret, Slater{SP}(state))
            println("HERE")
        end
    end

    ret
end

end # module States
