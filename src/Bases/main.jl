@reexport module Bases
export RefStates, AbstractState, Basis, SPBasis, RefState, MBBasis, Bra, overlap, ZeroState
# Basis interface
export basis, index, indexbasis, dim #, Base.indices

using Reexport: @reexport

abstract type AbstractState end
abstract type Basis <: AbstractState end
abstract type SPBasis <: Basis end
abstract type MBBasis <: Basis end
struct ZeroState <: AbstractState end

Base.:(==)(x::Basis, y::Basis) = false
Base.:(==)(x::B, y::B) where B<:Basis = index(x) == index(y)

Base.getindex(::Type{B}, ixs...) where B<:Basis = indexbasis(B, ixs...)
Base.getindex(::Type{B}, ixs::Array) where B<:Basis = map(i -> indexbasis(B, i), ixs)
Base.getindex(::Type{B}, r::Range{Int}) where B<:Basis = map(i -> indexbasis(B, i), r)

Base.start(::Type{B}) where B<:Basis = start(basis(B))
Base.next(::Type{B}, st) where B<:Basis = next(basis(B), st)
Base.done(::Type{B}, st) where B<:Basis = done(basis(B), st)

Base.endof(::Type{B}) where B<:Basis = dim(B)

Base.indices(::Type{B}) where B<:Basis = 1:dim(B)
basis(::Type{B}) where B<:Basis = map(x -> B[x], indices(B))

struct Bra{S<:AbstractState}; state::S end

Base.:(==)(b1::Bra{SP}, b2::Bra{SP}) where SP = b1.state == b2.state
Base.ctranspose(s::AbstractState) = Bra(s)
Base.:*(bra::Bra, ket::AbstractState) = overlap(bra.state, ket)
Base.:*(x::Number, a::ZeroState) = ZeroState()
Base.:*(a::ZeroState, x::Number) = ZeroState()

overlap(a::B, b::B) where B<:Basis = a == b
overlap(a::ZeroState, b::S) where S<:AbstractState = 0
overlap(a::S, b::ZeroState) where S<:AbstractState = 0
overlap(a::ZeroState, b::ZeroState) = 0

include("indexbasis.jl")
include("pairing.jl")
include("refstates.jl")
include("mbbasis.jl")
include("subbasis.jl")
include("states.jl")

function Base.show(io::IO, x::Union{PartHole{R}, <:SubBasis{PartHole{R}}}) where {SP, R<:RefState{SP}}
    x = convert(PartHole{R}, x)
    print(io, "Parts[", [string(indexp(R, i))*" " for (i, p) in enumerate(x.parts) if p]..., "]")
    print(" ")
    print(io, "Holes[", [string(indexh(R, i))*" " for (i, h) in enumerate(x.holes) if h]..., "]")
end

end # module States
