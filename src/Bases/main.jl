@reexport module Bases
#export AbstractBasis, ConcreteBasis, LeafBasis, RefStates, RefState

#=Modules=#  export RefStates
#=Types=#    export AbstractBasis, TensorBasis, Basis, RefState, Vacuum, Fermi
#=Shape=#    export rank, dim, fulldims
#=Indexing=# export index, linearindex, indexbasis
#=Algebra=#  export norm, overlap
#=Slater=#   export create, create!, annihil, annihil!
#=MBBasis=#  export occ,   unocc, nocc,   nunocc, isocc,  isunocc,
                    holes, parts, nholes, nparts, ishole, ispart,
                    spbasis

abstract type AbstractBasis
abstract type TensorBasis{Rank} <: AbstractBasis end
const Basis = TensorBasis{1}
#abstract type FockBasis{SPBasis<:TensorBasis} <: AbstractBasis end

rank(::Type{<:TensorBasis{N}}) where N = N
rank(b::TensorBasis) = rank(typeof(b))
dim(B::Type{<:TensorBasis}) = prod(fulldims(B))
# Good way to do this?
#fulldims(B::Type{<:Basis}) = (dim(B),)

Base.:(==)(x::B, y::B) where B<:AbstractBasis = index(x) == index(y)

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

norm(b::AbstractBasis) = one(Int)
norm(T, b::AbstractBasis) = one(T)
overlap(T, a::B, b::B) where B<:AbstractBasis = a == b ? norm(T, a) : zero(T)
overlap(a::B, b::B) where B<:AbstractBasis = overlap(Int, a, b)

include("indexing.jl")
include("iter.jl")
include("pairing.jl")
include("subbasis.jl")
include("product.jl")
include("mbbasis.jl")
include("slater.jl")

@defSub(Paired{P, L} <: Slater{Pairing{L}}) do s
    SP = Pairing{L}

    P == count(findall(s.occ)) do I
        s.occ[flipspin(SP[I])]
    end
end

end # module States
