@reexport module RefStates
export particles, holes

using ..RefState, ..SPState, ..level, ..iter

struct Vacuum{SP<:SPState} <: RefState end
struct Fermi{F, SP<:SPState} <: RefState end

particles(T::Type) = T()
particles(::Vacuum{SP}) where SP = iter(SP)
particles(::Fermi{F, SP}) where {F, SP} = (p for p in iter(SP) if level(p) <= F)

holes(T::Type) = T()
holes(::Vacuum) = ()
holes(::Fermi{F, SP}) where {F, SP} = (p for p in iter(SP) if level(p) > F)

end # module RefStates
