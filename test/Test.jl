module Test
using Base.Test
include("../src/ManyBody.jl")
using .ManyBody

include("basisindex.jl")
include("normordtest.jl")

end
