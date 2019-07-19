abstract type MBBasis{SPBasis<:AbstractBasis} <: AbstractBasis end
abstract type RefState{SPBasis<:AbstractBasis} end

spbasis(::Type{<:RefState{SPB}}) where SPB<:AbstractBasis = SPB
spbasis(::Type{<:MBBasis{SPB}}) where SPB<:AbstractBasis = SPB
spbasis(r::RefState) = spbasis(typeof(r))
spbasis(r::MBBasis) = spbasis(typeof(r))

# Must define
function isocc end

eachocc(x) = (p for p in elems(spbasis(x)) if isocc(x, p))
eachunocc(x) = (p for p in elems(spbasis(x)) if isunocc(x, p))
eachocc_index(x) = (SPB=spbasis(x); (i for i in eachindex(SPB) if isocc(x, (+SPB)[i])))
eachunocc_index(x) = (SPB=spbasis(x); (i for i in eachindex(SPB) if isunocc(x, (+SPB)[i])))
nocc(x) = length(eachocc(x))
nunocc(x) = length(eachunocc(x))
isunocc(x, p) = ~isocc(x, p)

const eachhole, eachpart = eachocc, eachunocc
const eachhole_index, eachpart_index = eachocc_index, eachunocc_index
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
eachocc(ref::RefStates.Vacuum) = ()
eachunocc(ref::RefStates.Vacuum) = elems(spbasis(ref))
nocc(::RefStates.Vacuum) = 0
nunocc(ref::RefStates.Vacuum) = dim(spbasis(ref))
isocc(ref::RefStates.Vacuum, s) = false
