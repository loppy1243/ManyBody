module Test

using Base.Test
include("../src/ManyBody.jl")
using .ManyBody

include("basisindex.jl")
include("normordtest.jl")
include("rlapply.jl")

function runall()
    basisindex()
    normordtest()
    rlapply()
end

end
