abstract type MBBasis{SPBasis<:AbstractBasis} <: AbstractBasis end
abstract type RefState{SPBasis<:AbstractBasis} end

spbasis(::Type{<:RefState{SPB}}) where SPB<:AbstractBasis = SPB
spbasis(::Type{<:MBBasis{SPB}}) where SPB<:AbstractBasis = SPB
spbasis(r::RefState) = spbasis(typeof(r))
spbasis(r::MBBasis) = spbasis(typeof(r))

# Must define
isocc(x) = MethodError(isocc, (x,)) |> throw
occ(x) = [p for p in spbasis(mb) if isocc(x, p)]
unocc(x) = [p for p in spbasis(mb) if isunocc(x, p)]
nocc(x) = length(occ(x))
nunocc(x) = length(unocc(x))
isunocc(x, p) = ~isocc(x, p)

const holes, parts = occ, unocc
const nholes, nparts = nocc, nunocc
const ishole, ispart = isocc, isunocc
Base.in(p, ref::RefState) = isocc(mb, p)
Base.in(p, mb::MBBasis) = isocc(mb, p)

module RefStates
    export Vacuum, Fermi

    using ..Bases
    
    struct Vacuum{SPB<:AbstractBasis} <: RefState{SPB} end
    struct Fermi{SPB<:AbstractBasis} <: RefState{SPB}; fermilevel::Int end

    Vacuum(SPB::Type{<:AbstractBasis}) = Vacuum{SPB}()
    Fermi{SPB}(F::Int) where SPB<:AbstractBasis = Fermi{SPB}(F)
end # module RefStates
occ(ref::RefStates.Vacuum) = Vector{spbasis(ref)}()
unocc(ref::RefStates.Vacuum) = collect(spbasis(ref))
nocc(::RefStates.Vacuum) = 0
nunocc(ref::RefStates.Vacuum) = dim(spbasis(ref))
isocc(ref::RefStates.Vacuum, s) = false
