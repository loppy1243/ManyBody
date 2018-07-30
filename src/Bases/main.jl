@reexport module Bases
export RefStates, State, Basis, SPBasis, RefState, MBBasis, Bra, overlap
# Basis interface
export index, indexbasis, dim #, Base.indices

using Reexport: @reexport

abstract type State end
abstract type Basis <: State end
abstract type SPBasis <: Basis end
abstract type MBBasis <: Basis end

Base.:(==)(x::Basis, y::Basis) = false
Base.:(==)(x::B, y::B) where B<:Basis = index(x) == index(y)

Base.getindex(::Type{B}, ixs...) where B<:Basis = indexbasis(B, ixs...)
Base.getindex(::Type{B}, ixs::Array) where B<:Basis = map(i -> indexbasis(B, i), ixs)

Base.start(::Type{B}) where B<:Basis = start(basis(B))
Base.next(::Type{B}, st) where B<:Basis = next(basis(B), st)
Base.done(::Type{B}, st) where B<:Basis = done(basis(B), st)

Base.indices(::Type{B}) where B<:Basis = 1:dim(B)

struct Bra{S<:State}; state::S end

Base.:(==)(::Bra, ::Bra) = false
Base.:(==)(b1::Bra{SP}, b2::Bra{SP}) where SP = b1.state == b2.state
Base.ctranspose(s::State) = Bra(s)
Base.:*(bra::Bra, ket::State) = overlap(bra.state, ket.state)

overlap(a, b) = MethodError(overlap, (a, b)) |> throw

include("indexbasis.jl")
include("pairing.jl")
include("refstates.jl")
include("mbbasis.jl")
include("subbasis.jl")

end # module States
