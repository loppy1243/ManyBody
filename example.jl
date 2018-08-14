module Exec

include("src/ManyBody.jl")

using Combinatorics: combinations
using JuliaUtil: cartesian_pow
using .ManyBody

Bases.@defSub NPairing{L, P} <: Bases.Slater{Bases.Pairing{L}} begin
    SP = Bases.Pairing{L}
    MB = Bases.Slater{SP}

    [MB(sps) for sps in combinations(SP, P)]
end

const LEVEL_SPACING = 1
const SPBASIS = Bases.Pairing{4}
const REFSTATE = RefStates.Fermi{2, SPBASIS}
#const MBBASIS = Bases.Paired{2, 4}
const MBBASIS = NPairing{4, 4}

f(g) = F64ActionOperator{1, MBBASIS}() do X
    sum(SPBASIS) do p
        LEVEL_SPACING*(level(p)-1)*A(p', p)(X)
    end
end

V(g) = F64ActionOperator{1, MBBASIS}() do X
    -g/2 * sum(cartesian_pow(1:nlevels(SPBASIS), Val{2})) do ls
        p, q, r, s = SPBASIS.((ls[1], ls[1], ls[2], ls[2]),
                              (SPINUP, SPINDOWN, SPINUP, SPINDOWN))
        A(p', q', s, r)(X)
    end
end

H(g) = f(g) + V(g)
main() = tabulate(H(1.0))

end # module Exec
