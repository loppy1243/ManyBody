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

Base.start(::Type{B}) where B<:AbstractBasis = start(basis(B))
Base.next(::Type{B}, st) where B<:AbstractBasis = next(basis(B), st)
Base.done(::Type{B}, st) where B<:AbstractBasis = done(basis(B), st)

Base.endof(::Type{B}) where B<:AbstractBasis = dim(B)

Base.indices(::Type{B}) where B<:AbstractBasis = 1:dim(B)
basis(::Type{B}) where B<:AbstractBasis = map(x -> B[x], indices(B))

include("indexbasis.jl")
include("subbasis.jl")
include("pairing.jl")
include("refstates.jl")
include("mbbasis.jl")

@defSubBasis Paired{F, L} <: Slater{Pairing{L}} begin
    SP = Pairing{L}

    ret = [Slater{SP}()]
    for ph = 1:L-F
        for lhs in combinations(1:F, ph), lps in combinations(F+1:L, ph)
            state = map(lhs, lps) do lh, lp
                [SP(lp, SPINUP),   SP(lh, SPINUP),
                 SP(lp, SPINDOWN), SP(lh, SPINDOWN),
                 SP(lp, SPINUP),   SP(lh, SPINDOWN),
                 SP(lp, SPINDOWN), SP(lh, SPINUP)  ]
            end |> x -> reduce(vcat, x)

            push!(ret, Slater{SP}(state))
        end
    end

    ret
end

end # module States
