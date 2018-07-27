@reexport module RefStates
export particles, holes

using ..RefState, ..SPState, ..level, ..states

struct Vacuum{SP<:SPState} <: RefState end
struct Fermi{F, SP<:SPState} <: RefState end

parts(::Type{<:Vacuum}) = states(SP)
parts(::Type{Fermi{F, SP}}) where {F, SP} = (p for p in states(SP) if level(p) <= F)

holes(::Type{<:Vacuum}) = ()
holes(::Type{Fermi{F, SP}}) where {F, SP} = (p for p in states(SP) if level(p) > F)

nparts(::Type{R}) where R<:RefState = length(parts(R))
nparts(::Type{Vacuum{SP}}) where SP = nstates(SP)
nparts(::Type{Fermi{F, PairingParticle{L}}}) where {F, L} = 2(L-F)

nholes(::Type{R}) where R<:RefState = length(holes(r))
nholes(::Type{<:Vacuum}) = 0
nholes(::Type{Fermi{F, PairingParticle{L}}}) where {F, L} = 2F

pnum(::Type{R}, s) where R<:RefState = findin(parts(s), s)
pnum(::Type{Vacuum{SP}}, s::SP) where SP = snum(s)
pnum(::Type{Fermi{F, PairingParticle{L}}}, s::PairingParticle{L}) where {F, L} =
    pnum(s) - nholes(s) + 1

hnum(::Type{R}, s) where R<:RefState = findin(holes(s), s)
# Don't know if this should even be defined
hnum(r::Type{Vacuum{SP}}, s::SP) where SP = 0
hnum(::Type{Fermi{F, PairingParticle{L}}}, s::PairingParticle{L}) where {F, L} =
    pnum(s)

isocc(::Type{R}, s) where {R<:RefState, SP} = s in holes(R)
isocc(::Type{Vacuum{SP}}, s::SP) where SP = false

end # module RefStates
