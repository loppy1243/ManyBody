@reexport module RefStates
export RefState, Vacuum, Fermi, parts, holes, nparts, pindex, hindex, occ, unocc, n_occ,
       n_unocc, isocc, ishole, ispart, nholes, nparts, indexh, indexp, fermilevel

using ..Bases

abstract type RefState{SP<:SPBasis} end

struct Vacuum{SP<:SPBasis} <: RefState{SP} end
struct Fermi{F, SP<:SPBasis} <: RefState{SP} end

fermilevel(::Type{<:Fermi{F}}) where F = F

parts(::Type{<:Vacuum}) = states(SP)
parts(::Type{Fermi{F, SP}}) where {F, SP} = [p for p in SP if level(p) > F]
const unocc = parts

holes(::Type{<:Vacuum}) = ()
holes(::Type{Fermi{F, SP}}) where {F, SP} = [p for p in SP if level(p) <= F]
const occ = holes

nparts(::Type{R}) where R<:RefState = length(parts(R))
nparts(::Type{Vacuum{SP}}) where SP = nstates(SP)
nparts(::Type{Fermi{F, Pairing{L}}}) where {F, L} = 2(L-F)
const n_unocc = nparts

nholes(::Type{R}) where R<:RefState = length(holes(R))
nholes(::Type{<:Vacuum}) = 0
nholes(::Type{Fermi{F, <:Pairing}}) where F = 2F
const n_occ = nholes

pindex(::Type{R}, s) where R<:RefState = findfirst(parts(R), s)
pindex(::Type{Vacuum{SP}}, s::SP) where SP = index(s)
pindex(::Type{R}, s::SP) where {SP<:Pairing, R<:Fermi{<:Any, SP}} = index(s) - nholes(R)

hindex(::Type{R}, s::SP) where {SP, R<:RefState{SP}} = findfirst(holes(R), s)
# Don't know if this should even be defined
hindex(::Type{Vacuum{SP}}, s::SP) where SP = 0
hindex(::Type{R}, s::SP) where {SP<:Pairing, R<:Fermi{<:Any, SP}} = index(s)

indexp(::Type{R}, n::Int) where R<:RefState = parts(R)[n]
indexp(::Type{Vacuum{SP}}, n::Int) where SP = SP[n]
indexp(::Type{Fermi{F, Pairing{L}}}, n::Int) where {F, L} =
    Pairing{L}[nholes(Fermi{F, Pairing{L}}) + n]

indexh(::Type{R}, n::Int) where R<:RefState = holes(SP)[n]
indexh(::Type{<:Vacuum}, n::Int) = 0
indexh(::Type{Fermi{F, SP}}, n::Int) where {F, SP<:Pairing} = SP[n]

isocc(::Type{R}, s::SP) where {SP, R<:RefState{SP}} = s in holes(R)
isocc(::Type{Vacuum{SP}}, s::SP) where SP = false
isocc(::Type{R}, s::SP) where {F, SP<:Pairing, R<:Fermi{F, SP}} =
    level(s) <= F
const ishole = isocc

ispart(R, x) = ~ishole(R, x)

end # module RefStates
