@reexport module Bases
export RefStates, State, Basis, SPBasis, RefState, Index, MBBasis, Bra, overlap, indexes, dim,
       indexbasis, index

using Reexport: @reexport

abstract type State end
abstract type RefState end
abstract type Basis <: State end
abstract type SPBasis <: Basis end

Base.getindex(::Type{B}, ixs) where B<:Basis = indexbasis(B, ixs)
Base.getindex(::Type{B}, ixs...) where B<:Basis = indexbasis(B, ixs)

Base.start(::Type{B}) where B<:Basis = start(basis(B))
Base.next(::Type{B}, st) where B<:Basis = next(basis(B), st)
Base.done(::Type{B}, st) where B<:Basis = done(basis(B), st)

struct Index{SP<:SPBasis} <: SPBasis
    index::Int
end
Index(s::SPBasis) = Index{typeof(s)}(index(s))

Base.indices(::Type{B}) where B<:Basis = 1:dim(B)
indexes(::Type{B}) where B<:Basis = map(Index{B}, indices(B))

index(s::Index) = s.index
indexbasis(::Type{S}, s::Int) where S<:Index = S(s)
basis(::Type{Index{B}}) where B = indexes(B)
dim(::Type{Index{B}}) where B = dim(B)

Base.convert(::Type{Index{B}}, s::B) where B = Index{B}(index(s))
Base.convert(::Type{Index}, s::Basis) = Index{typeof(s)}(index(s))
Base.convert(::Type{B}, s::Index{B}) where B<:Basis = B(s)

SP(s::Index{SP}) where SP<:SPBasis = SP[s.index]

struct Bra{S<:State}; state::S end

Base.:(==)(::Bra, ::Bra) = false
Base.:(==)(b1::Bra{SP}, b2::Bra{SP}) where SP = b1.state == b2.state
Base.ctranspose(s::State) = Bra(s)
Base.:*(bra::Bra, ket::State) = overlap(bra.state, ket.state)

overlap(a, b) = MethodError(overlap, (a, b)) |> throw

include("pairing.jl")
include("refstates.jl")
include("mbbasis.jl")

end # module States
