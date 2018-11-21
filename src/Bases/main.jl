@reexport module Bases
#export AbstractBasis, ConcreteBasis, LeafBasis, RefStates, RefState

#=Modules=#  export RefStates
#=Types=#    export AbstractBasis, TensorBasis, Basis, RefState, Vacuum, Fermi
#=Shape=#    export rank, dim, fulldims
#=Indexing=# export indexer, linearindexer, index, linearindex, indexbasis
#=Algebra=#  export norm, overlap
#=Sub=#      export superbasis, supbasis, superindex, supindex, superelem, supelem, subindexmap
#=Slater=#   export create, create!, annihil, annihil!, createsgn, annihilsgn
#=MBBasis=#  export occ,   unocc, occinds,  unoccinds, nocc,   nunocc, isocc,  isunocc,
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
    ret[linearindex(b)] = oneunit(T)

    ret
end
Base.convert(::Type{Array{T}}, b::TensorBasis) where T = convert(Array{T, rank(b)}, b)
(A::Type{<:Array})(b::AbstractBasis) = convert(A, b)

norm(b::AbstractBasis) = 1
overlap(a::B, b::B) where B<:AbstractBasis = a == b ? norm(a) : 0

include("mbbasis.jl")
include("subbasis.jl")

include("pairing.jl")
include("product.jl")
include("slater.jl")

include("indexing.jl")
include("iter.jl")

@defSub Paired{A, L} <: Slater{Pairing{L}, 2} begin s
    SPI = indexer(Pairing{L})

    cnt = 0
    for I in occinds(s)
        cnt += flipspin(SPI[I]) in s || return false
    end
    cnt == A
end

@defSub MBPairing{A, L} <: Slater{Pairing{L}, 2} begin s
    nocc(s) == A
end

end # module States
