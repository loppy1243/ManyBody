@reexport module States
export RefStates, State, SPState, RefState, MBState, pnum, nlevels, level, spin, iter

using Reexport: @reexport

abstract type State
abstract type SPState <: State
abstract type RefState <: State

include("mbstate.jl")
include("pairing.jl")
include("refstates.jl")

end # module States
