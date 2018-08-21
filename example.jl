module Exec

#include("src/ManyBody.jl")

using ManyBody
using ManyBody.Operators: A
using Combinatorics: combinations
using JuliaUtil: cartesian_pow

Bases.@defSub NPairing{L, P} <: Bases.Slater{Bases.Pairing{L}} begin
    SP = Bases.Pairing{L}
    MB = Bases.Slater{SP}

    [MB(sps) for sps in combinations(SP, P)]
end
Bases.IndexType(::Type{<:NPairing}) = Bases.IndexTypes.Linear()

const LEVEL_SPACING = 1
const SPBASIS = Bases.Pairing{4}
const REFSTATE = RefStates.Fermi{2, SPBASIS}
const MBBASIS = Bases.Paired{2, 4}
#const MBBASIS = NPairing{4, 4}

f(g) = F64FunctionOperator{MBBASIS}() do X, Y
    sum(SPBASIS) do p
        LEVEL_SPACING*(level(p)-1)*(X'A(p', p)(Y))
    end
end

V(g) = F64FunctionOperator{MBBASIS}() do X, Y
    -g/2 * sum(cartesian_pow(1:nlevels(SPBASIS), Val{2})) do ls
        p, q, r, s = SPBASIS.((ls[1], ls[1], ls[2], ls[2]),
                              (SPINUP, SPINDOWN, SPINUP, SPINDOWN))
        X'A(p', q', s, r)(Y)
    end
end

H(g) = f(g) + V(g)
main() = tabulate(H(1.0))

end # module Exec
main()
