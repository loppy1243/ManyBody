module Exec

using ManyBody
using ManyBody.Operators: @A

Bases.@defSub NPairing{L, P} <: Bases.Slater{Bases.Pairing{L}} do b
    ps = occ(b)

    length(ps) && all(p -> p.level <= L, occ(b))
end

const LEVEL_SPACING = 1
const SPBASIS = Bases.Pairing{4}
const REFSTATE = RefStates.Fermi{SPBASIS, 2}
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
Exec.main()
