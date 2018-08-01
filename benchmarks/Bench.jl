module Bench
using BenchmarkTools
include("../src/ManyBody.jl")
using .ManyBody

include("rlexpect.jl")

end # module Bench
