module Test

using Base.Test
include("../src/ManyBody.jl")
using .ManyBody

include("basisindex.jl")
include("normordtest.jl")

function runall()
    basisindex()
    normordtest()
end

end
