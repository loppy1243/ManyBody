@reexport module States
export RefStates, State, SPState, RefState, MBState, snum, nlevels, level, spin, iter, Bra,
       overlap

using Reexport: @reexport

abstract type State
abstract type SPState <: State
abstract type RefState <: State

struct IntState{SP<:SPState} <: SPState
    num::Int
end
IntState(s::SPState) = IntState{typeof(s)}(s)
snum(IntState) = s.num
states(::Type{IntState{SP}}) where SP = (IntState{SP}(n) for n = 1:nstates(SP))
nstates(::Type{IntState{SP}}) where SP = nstates(SP)

Base.convert(::Type{IntState{SP}}, s::SP) where SP = IntState{SP}(snum(s))
Base.convert(::Type{IntState}, s::SPState) = IntState{typeof(s)}(snum(s))
Base.convert(::Type{SP}, s::IntState{SP}) where SP = SP(s)
Base.convert(::Type{SPState}, s::IntState{SP}) where SP = SP(s)

SP(s::IntState{SP}) where SP<:SPState = particle(SP, s.num)

struct Bra{S<:State}; state::S end
Bra(s::State) = Bra{typeof(s)}(s)

Base.:(==)(::Bra, ::Bra) = false
Base.:(==)(b1::Bra{SP}, b2::Bra{SP}) where SP = b1 == b2
Base.ctranspose(s::State) = Bra(s)
Base.:*(bra::Bra, ket::State) = overlap(bra.state, ket.state)

overlap(a, b) = MethodError(overlap, (a, b)) |> throw

include("mbstate.jl")
include("pairing.jl")
include("refstates.jl")

end # module States
