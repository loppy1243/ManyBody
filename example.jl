module Exec

include("src/ManyBody.jl")

using JuliaUtil: cartesian_pow
using .ManyBody

const LEVEL_SPACING = 1
const SPBASIS = Bases.Pairing{4}
const REFSTATE = RefStates.Fermi{2, SPBASIS}
const MBBASIS = Bases.Paired{2, 4}
#const INDICES = indices(SPBASIS)
#const PARTS = parts(REFSTATE)
#const HOLES = holes(REFSTATE)

f(g) = (X, Y) -> sum(SPBASIS) do p
    sgn, Y2 = A(p', p)(Y)
    (LEVEL_SPACING*(level(p)-1))*sgn*(X'Y2)
end

V(g) = (X, Y) -> sum(cartesian_pow(1:nlevels(SPBASIS), Val{2})) do ls
    p, q, r, s = SPBASIS.((ls[1], ls[1], ls[2], ls[2]),
                          (SPINUP, SPINDOWN, SPINUP, SPINDOWN))

    sgn, Y2 = A(p', q', s, r)(Y)
    -g/2*sgn*(X'Y2)
end

H(g) = Operator{2, SPBASIS}() do X, Y
    f(g)(X, Y) + V(g)(X, Y)
end

main() = tabulate(H(1.0))

end # module Exec
