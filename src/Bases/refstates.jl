export RefStates, RefState, Vacuum, Fermi, parts, holes, nparts, pindex, hindex, occ, unocc,
       n_occ, n_unocc, isocc, isunocc, ishole, ispart, nholes, nparts
using ..ManyBody

abstract type RefState{SPB<:AbstractBasis} end

module RefStates
    export Vacuum, Fermi

    using ..Bases
    
    struct Vacuum{SPB<:AbstractBasis} <: RefState{SPB} end
    struct Fermi{SPB<:AbstractBasis} <: RefState{SPB}; fermilevel::Int end

    Vacuum(SPB::Type{<:AbstractBasis}) = Vacuum{SPB}()
    Fermi{SPB}(F::Int) where SPB<:AbstractBasis = Fermi{SPB}(F)
end # module RefStates

spbasis(::Type{<:RefState{SPB}}) where SPB<:AbstractBasis = SPB
spbasis(r::RefState) = spbasis(typeof(r))

unocc(ref::RefStates.Vacuum) = collect(spbasis(ref))
unocc(ref::RefStates.Fermi) = [p for p in spbasis(ref) if p.level > ref.fermilevel]
const parts = unocc

occ(ref::RefStates.Vacuum) = Vector{spbasis(ref)}()
occ(ref::RefStates.Fermi) = [p for p in spbasis(ref) if p.level <= ref.fermilevel]
const holes = occ

nunocc(ref::RefState) = length(parts(ref))
nunocc(ref::RefStates.Vacuum) = dim(spbasis(ref))
const nparts = nunocc

nocc(::RefState) = length(holes(R))
nocc(::RefStates.Vacuum) = 0
const nholes = nocc

isocc(ref::RefState{SPB}, b::SPB) where SPB<:AbstractBasis = b in holes(ref)
isocc(ref::RefState, b::AbstractBasis) = isocc(ref, convert(spbasis(ref), b))
isocc(ref::RefStates.Vacuum, s) = false
const ishole = isocc

isunocc(r, s) = ~isocc(r, s)
const ispart = isunocc
