export RefStates, RefState, Vacuum, Fermi, parts, holes, nparts, pindex, hindex, occ, unocc,
       n_occ, n_unocc, isocc, ishole, ispart, nholes, nparts, fermilevel

abstract type RefState{B<:ConcreteBasis} end

module RefStates
    export Vacuum, Fermi

    using ..Bases
    
    struct Vacuum{B<:ConcreteBasis} <: RefState{B} end
    struct Fermi{B<:ConcreteBasis, F} <: RefState{B} end
end # module RefStates

ManyBody.basistype(::Type{<:RefState{B}}) where B<:ConcreteBasis = B

fermilevel(::RefStates.Fermi{<:Any, F}) where F = F

parts(v::RefStates.Vacuum) = basis(basistype(v))
parts(f::RefStates.Fermi) = [p for p in basistype(f) if level(p) > fermilevel(f)]
const unocc = parts

#holes(v::RefStates.Vacuum) = basistype(v)[]
holes(::RefStates.Vacuum) = ()
holes(f::RefStates.Fermi) = [p for p in SP if level(p) <= fermilevel(f)]
const occ = holes

nparts(r::RefState) = length(parts(r))
nparts(v::RefStates.Vacuum) = dim(basistype(v))
nparts(f::RefStates.Fermi{<:Pairing}) = 2(nlevels(basistype(f)) - fermilevel(f))
const n_unocc = nparts

nholes(::RefState) = length(holes(R))
nholes(::RefStates.Vacuum) = 0
nholes(f::RefStates.Fermi{<:Pairing}) = 2fermilevel(f)
const n_occ = nholes

isocc(r::RefState, s::Wrapped) = isocc(r, inner(s))

isocc(r::RefState{B}, s::B) where B<:ConcreteBasis = s in holes(r)
isocc(v::RefStates.Vacuum, s) = false
isocc(f::RefStates.Fermi{B}, s::B) where B<:Pairing =
    level(s) <= fermilevel(f)
isocc(f::RefStates.Fermi{B}, s::MaybeIndex{B}) where B<:Pairing = index(s) <= 2fermilevel(f)
const ishole = isocc

ispart(R, x) = ~ishole(R, x)
