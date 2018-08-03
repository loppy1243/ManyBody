export RefStates, RefState, Vacuum, Fermi, parts, holes, nparts, pindex, hindex, occ, unocc,
       n_occ, n_unocc, isocc, ishole, ispart, nholes, nparts, fermilevel

abstract type RefState{B<:AbstractBasis} end

module RefStates
    export Vacuum, Fermi

    using ..Bases
    
    struct Vacuum{B<:AbstractBasis} <: RefState{B} end
    struct Fermi{F, B<:AbstractBasis} <: RefState{B} end
end # module RefStates

fermilevel(::Type{<:RefStates.Fermi{F}}) where F = F

parts(::Type{<:RefStates.Vacuum}) = states(SP)
parts(::Type{RefStates.Fermi{F, SP}}) where {F, SP} = [p for p in SP if level(p) > F]
const unocc = parts

holes(::Type{<:RefStates.Vacuum}) = ()
holes(::Type{RefStates.Fermi{F, SP}}) where {F, SP} = [p for p in SP if level(p) <= F]
const occ = holes

nparts(::Type{R}) where R<:RefState = length(parts(R))
nparts(::Type{RefStates.Vacuum{SP}}) where SP = nstates(SP)
nparts(::Type{RefStates.Fermi{F, Pairing{L}}}) where {F, L} = 2(L-F)
const n_unocc = nparts

nholes(::Type{R}) where R<:RefState = length(holes(R))
nholes(::Type{<:RefStates.Vacuum}) = 0
nholes(::Type{RefStates.Fermi{F, <:Pairing}}) where F = 2F
const n_occ = nholes

#pindex(::Type{R}, s) where R<:RefState = findfirst(parts(R), s)
#pindex(::Type{RefStates.Vacuum{SP}}, s::SP) where SP = index(s)
#pindex(::Type{R}, s::SP) where {SP<:Pairing, R<:RefStates.Fermi{<:Any, SP}} = index(s) - nholes(R)
#
#hindex(::Type{R}, s::SP) where {SP, R<:RefState{SP}} = findfirst(holes(R), s)
## Don't know if this should even be defined
#hindex(::Type{RefStates.Vacuum{SP}}, s::SP) where SP = 0
#hindex(::Type{R}, s::SP) where {SP<:Pairing, R<:RefStates.Fermi{<:Any, SP}} = index(s)
#
#indexp(::Type{R}, n::Int) where R<:RefState = parts(R)[n]
#indexp(::Type{RefStates.Vacuum{SP}}, n::Int) where SP = SP[n]
#indexp(::Type{RefStates.Fermi{F, Pairing{L}}}, n::Int) where {F, L} =
#    Pairing{L}[nholes(RefStates.Fermi{F, Pairing{L}}) + n]
#
#indexh(::Type{R}, n::Int) where R<:RefState = holes(SP)[n]
#indexh(::Type{<:RefStates.Vacuum}, n::Int) = 0
#indexh(::Type{RefStates.Fermi{F, SP}}, n::Int) where {F, SP<:Pairing} = SP[n]

isocc(::Type{R}, s::SP) where {SP, R<:RefState{SP}} = s in holes(R)
isocc(::Type{RefStates.Vacuum{SP}}, s::SP) where SP = false
isocc(::Type{R}, s::SP) where {F, SP<:Pairing, R<:RefStates.Fermi{F, SP}} =
    level(s) <= F
const ishole = isocc

ispart(R, x) = ~ishole(R, x)
