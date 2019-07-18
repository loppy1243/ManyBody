@reexport module Bases
#export AbstractBasis, ConcreteBasis, LeafBasis, RefStates, RefState

#=Modules=#  export RefStates
#=Types=#    export AbstractBasis, TensorBasis, Basis, MBBasis, RefState, Vacuum, Fermi
#=Shape=#         export rank, dim, fulldims
#=Indexing/Iter=# export index, indexbasis, elems
#=Algebra=#       export norm, overlap
#=Sub=#           export superbasis, supbasis, superindex, supindex, superelem, supelem,
                         subindexmap
#=Slater=#        export create, create!, annihil, annihil!, createsgn, annihilsgn
#=MBBasis=#       export occ,   unocc, occinds,  unoccinds, nocc,   nunocc, isocc,  isunocc,
                         holes, parts, holeinds, partinds,  nholes, nparts, ishole, ispart,
                         spbasis

abstract type AbstractBasis end
abstract type TensorBasis{Rank} <: AbstractBasis end
const Basis = TensorBasis{1}
#abstract type FockBasis{SPBasis<:TensorBasis} <: AbstractBasis end

rank(::Type{<:TensorBasis{N}}) where N = N
rank(b::TensorBasis) = rank(typeof(b))
dim(B::Type{<:TensorBasis}) = prod(fulldims(B))
# Good way to do this?
#fulldims(B::Type{<:Basis}) = (dim(B),)

Base.:(==)(x::B, y::B) where B<:AbstractBasis = index(x) == index(y)
Base.in(::B, ::Type{B}) where B<:AbstractBasis = true
Base.in(::AbstractBasis, ::Type{<:AbstractBasis}) = false

function Base.convert(::Type{Array{T, N}}, b::TensorBasis{N}) where {T, N}
    ret = zeros(T, fulldims(typeof(b)))
    ret[b] = oneunit(T)

    ret
end
function Base.convert(::Type{Vector{T}}, b::AbstractBasis) where T
    ret = zeros(T, dim(typeof(b)))
    ret[-b] = oneunit(T)

    ret
end
Base.convert(::Type{Array{T}}, b::TensorBasis) where T = convert(Array{T, rank(b)}, b)
(A::Type{<:Array})(b::AbstractBasis) = convert(A, b)

norm(b::AbstractBasis) = 1
overlap(a::B, b::B) where B<:AbstractBasis = a == b ? norm(a) : 0

include("arrayinterface.jl")

include("mbbasis.jl")
include("subbasis.jl")

include("pairing.jl")
include("product.jl")
include("slater.jl")

@defSub Paired{A, L} <: Slater{Pairing{L}, 2} begin s
    cnt = 0
    for p in occ(s)
        cnt += flipspin(p) in s || return false
    end
    cnt == A
end

@defSub MBPairing{A, L} <: Slater{Pairing{L}, 2} begin s
    nocc(s) == A
end

Base.show(io::IO, mime::MIME"text/plain", B::Type{<:AbstractBasis}) = show(io, mime, +B)

end # module States
