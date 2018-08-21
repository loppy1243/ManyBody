module ManyBody
export Bases, States

using Reexport: @reexport

abstract type AbstractState end

include("util.jl")
using .@disallow
include("SpinMod.jl")
include("Bases/main.jl")

basistype(x) = basistype(typeof(x))
@disallow basistype(::Type)
basistype(x::AbstractBasis) = basistype(typeof(x))
basistype(B::Type{<:AbstractBasis}) = B

include("States.jl")
include("Operators.jl")

end # module ManyBody
