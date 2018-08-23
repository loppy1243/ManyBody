module ManyBody
export Bases, States, Operators, basistype, reptype

using Reexport: @reexport

include("util.jl")

abstract type AbstractState end

basistype(x) = basistype(typeof(x))
@disallow basistype(::Type)

reptype(x) = reptype(typeof(x))
@disallow reptype(::Type)

include("SpinMod.jl")
include("Bases/main.jl")
include("States.jl")
include("Operators.jl")

end # module ManyBody
