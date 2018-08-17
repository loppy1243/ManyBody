module ManyBody
export Bases, States

using Reexport: @reexport

abstract type AbstractState end

include("SpinMod.jl")
include("Bases/main.jl")
include("States.jl")
include("Operators.jl")

end # module ManyBody
