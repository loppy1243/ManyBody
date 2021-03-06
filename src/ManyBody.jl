module ManyBody
export Bases, States, Operators, basistype, reptype

using Reexport: @reexport

abstract type AbstractState end

include("SpinMod.jl")
include("Bases/main.jl")
include("Operators.jl")
include("Hamiltonians.jl"); export Hamiltonians

end # module ManyBody
