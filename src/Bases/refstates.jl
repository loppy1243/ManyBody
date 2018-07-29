@reexport module RefStates
export Vacuum, Fermi, parts, holes, nparts, pindex, hindex, noccs, isocc, nholes, nparts,
       indexh, indexp

using ..Bases

struct Vacuum{SP<:SPBasis, N} <: RefState end
struct Fermi{F, SP<:SPBasis} <: RefState end

parts(::Type{<:Vacuum}) = states(SP)
parts(::Type{Fermi{F, SP}}) where {F, SP} = [p for p in SP if level(p) <= F]

holes(::Type{<:Vacuum}) = ()
holes(::Type{Fermi{F, SP}}) where {F, SP} = [p for p in SP if level(p) > F]

nparts(::Type{R}) where R<:RefState = length(parts(R))
nparts(::Type{<:Vacuum{SP}}) where SP = nstates(SP)
nparts(::Type{Fermi{F, Pairing{L}}}) where {F, L} = 2(L-F)

nholes(::Type{R}) where R<:RefState = length(holes(R))
nholes(::Type{<:Vacuum}) = 0
nholes(::Type{Fermi{F, <:Pairing}}) where F = 2F

pindex(::Type{R}, s) where R<:RefState = findin(parts(R), s)
pindex(::Type{<:Vacuum{SP}}, s::SP) where SP = index(s)
pindex(::Type{R}, s::SP) where {SP<:Pairing, R<:Fermi{<:Any, SP}} = index(s) - nholes(R)

hindex(::Type{R}, s) where R<:RefState = findin(holes(R), s)
# Don't know if this should even be defined
hindex(::Type{<:Vacuum{SP}}, s::SP) where SP = 0
hindex(::Type{R}, s::SP) where {SP<:Pairing, R<:Fermi{<:Any, SP}} = index(s)

indexp(::Type{<:Vacuum{SP}}, n::Int) where SP = SP[n]
indexp(::Type{R}, n::Int) where R<:Fermi = parts(R)[n]
indexp(::Type{Fermi{F, Pairing{L}}}, n::Int) where {F, L} =
    Pairing{L}[nholes(Fermi{F, Pairing{L}}) + n]

indexh(::Type{<:Vacuum}, n::Int) = 0
indexh(::Type{SP}, n::Int) where SP<:Fermi = holes(SP)[n]
indexh(::Type{Fermi{F, SP}}, n::Int) where {F, SP<:Pairing} = SP[n]

noccs(::Type{Vacuum{SP, N}}) where {SP, N} = N
noccs(::Type{<:Fermi{F}}) where F = 2F

isocc(::Type{R}, s) where R<:RefState = s in holes(R)
isocc(::Type{<:Vacuum{SP}}, s::SP) where SP = false
isocc(::Type{R}, s::SP) where {F, SP<:Pairing, R<:Fermi{F, SP}} =
    level(s) <= F

end # module RefStates
