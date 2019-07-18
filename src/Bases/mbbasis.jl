abstract type MBBasis{SPBasis<:AbstractBasis} <: AbstractBasis end
abstract type RefState{SPBasis<:AbstractBasis} end

spbasis(::Type{<:RefState{SPB}}) where SPB<:AbstractBasis = SPB
spbasis(::Type{<:MBBasis{SPB}}) where SPB<:AbstractBasis = SPB
spbasis(r::RefState) = spbasis(typeof(r))
spbasis(r::MBBasis) = spbasis(typeof(r))

# Must define
isocc(x) = MethodError(isocc, (x,)) |> throw

occ(x) = [p for p in elems(spbasis(x)) if isocc(x, p)]
occinds(x) = [SPB=spbasis(x); (i for i in eachindex(SPB) if isocc(x, (+SPB)[i]))]
unocc(x) = [p for p in elems(spbasis(x)) if isunocc(x, p)]
unoccinds(x) = [SPB=spbasis(x); (i for i in eachindex(SPB) if isunocc(x, (+SPB)[i]))]
nocc(x) = length(occ(x))
nunocc(x) = length(unocc(x))
isunocc(x, p) = ~isocc(x, p)

const holes, parts = occ, unocc
const holeinds, partinds = occinds, unoccinds
const nholes, nparts = nocc, nunocc
const ishole, ispart = isocc, isunocc
Base.in(p, ref::RefState) = isocc(mb, p)
Base.in(p, mb::MBBasis) = isocc(mb, p)

module RefStates
    export Vacuum, Fermi

    using ..Bases
    
    struct Vacuum{SPB<:AbstractBasis} <: RefState{SPB} end
    struct Fermi{SPB<:AbstractBasis} <: RefState{SPB}; fermilevel::Int end
end # module RefStates
occ(ref::RefStates.Vacuum) = Vector{spbasis(ref)}()
unocc(ref::RefStates.Vacuum) = elems(spbasis(ref))
nocc(::RefStates.Vacuum) = 0
nunocc(ref::RefStates.Vacuum) = dim(spbasis(ref))
isocc(ref::RefStates.Vacuum, s) = false
