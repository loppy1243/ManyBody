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

parts(ref::RefStates.Vacuum) = collect(spbasis(ref))
parts(ref::RefStates.Fermi) = [p for p in spbasis(ref) if p.level > ref.fermilevel]
const unocc = parts

holes(ref::RefStates.Vacuum) = Vector{spbasis(ref)}()
holes(ref::RefStates.Fermi) = [p for p in spbasis(ref) if p.level <= ref.fermilevel]
const occ = holes

nparts(ref::RefState) = length(parts(ref))
nparts(ref::RefStates.Vacuum) = dim(spbasis(ref))
const n_unocc = nparts

nholes(::RefState) = length(holes(R))
nholes(::RefStates.Vacuum) = 0
const n_occ = nholes

isocc(ref::RefState{SPB}, b::SPB) where SPB<:AbstractBasis = b in holes(ref)
isocc(ref::RefState, b::AbstractBasis) = isocc(ref, convert(spbasis(ref), b))
isocc(ref::RefStates.Vacuum, s) = false
const ishole = isocc

isunocc(r, s) = ~isocc(r, s)
const ispart = isunocc
