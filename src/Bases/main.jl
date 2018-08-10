@reexport module Bases
export RefStates, AbstractState, AbstractBasis, RefState, overlap
# Basis interface
export basis, index, indexbasis, dim

using ..AbstractState
using Combinatorics: combinations
using Reexport: @reexport

abstract type AbstractBasis <: AbstractState end

Base.:(==)(x::AbstractBasis, y::AbstractBasis) = false
Base.:(==)(x::B, y::B) where B<:AbstractBasis = index(x) == index(y)

Base.getindex(::Type{B}, ixs...) where B<:AbstractBasis = indexbasis(B, ixs...)
Base.getindex(::Type{B}, ixs::Array) where B<:AbstractBasis = map(i -> indexbasis(B, i), ixs)
Base.getindex(::Type{B}, r::AbstractRange{Int}) where B<:AbstractBasis = map(i -> indexbasis(B, i), r)
Base.firstindex(::Type{B}) where B<:AbstractBasis = 1
Base.lastindex(::Type{B}) where B<:AbstractBasis = dim(B)

Base.iterate(::Type{B}, i=1) where B<:AbstractBasis = i > dim(B) ? nothing : (B[i], i+1)
Base.IteratorSize(::Type{B}) where B<:AbstractBasis = Base.HasLength()
Base.IteratorEltype(::Type{B}) where B<:AbstractBasis = Base.HasEltype()

Base.length(::Type{B}) where B<:AbstractBasis = dim(B)
Base.eltype(::Type{B}) where B<:AbstractBasis = B

basis(::Type{B}) where B<:AbstractBasis = B[1:dim(B)]

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
            for (lh, lp) in zip(lhs, lps), s in SPINS
                annihil!(state, SP(lh, s))
                create!(state, SP(lp, s))
            end

            push!(ret, state)
        end
    end

    ret
end

Base.convert(::Type{B}, s::SubBasis{Index{B}}) = convert(B, s.state)
Base.convert(::Type{B}, s::Index{SB}) where SB<:SubBasis{B} = convert(B, SB[index(s)])
function Base.convert(::Type{V}, s::AbstractBasis) where V<:AbstractVector
    ret = V(undef, dim(typeof(s)))
    ret .= zero(eltype(V))
    ret[index(s)] = oneunit(eltype(V))

    ret
end
function Base.convert(::Type{B}, v::AbstractVector) where B<:AbstractBasis
    nzs = v .!= zero(eltype(v))
    i = findfirst(!iszero, nzs)
    if count(nzs) != 1 || v[i] != oneunit(eltype(v))
        InexactError() |> throw
    end

    B[i]
end


struct Bra{S}; state::S end
struct Zero end

Base.:*(a::Bra{B}, b::B) where B<:AbstractBasis = a == b
Base.:*(::Bra{Zero}, ::Union{<:AbstractVector, <:AbstractBasis}) = 0
Base.:*(::Union{<:Adjoint{<:Any, <:AbstractVector}, Bra{<:AbstractBasis}}, ::Zero) = 0
Base.:*(::Bra{Zero}, ::Zero) = 0
Base.:*(a::Bra{<:AbstractBasis}, b::AbstractVector) = b[index(a)]
Base.:*(a::Adjoint{<:Any, <:AbstractVector}, b::AbstractBasis) = conj(b[index(a)])

end # module States
