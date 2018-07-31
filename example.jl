@everywhere module Exec

include("src/ManyBody.jl")

using Loppy.Util: cartesian_pow
using .ManyBody

const LEVEL_SPACING = 1
const SPBASIS = Bases.Pairing{4}
const REFSTATE = RefStates.Fermi{2, SPBASIS}
const MBBASIS = Bases.PHPaired{2, 4}
#const INDICES = indices(SPBASIS)
#const PARTS = parts(REFSTATE)
#const HOLES = holes(REFSTATE)

const f_mb = @Operator(MBBASIS) do X, Y
    println("f: ", (index(X), index(Y)))
    sum(SPBASIS) do p
        (X == Y && p in X)*LEVEL_SPACING*(level(p)-1)
    end
end

V_sp(g) = @Operator(SPBASIS) do p, q, r, s
    if level(p) == level(q) && level(r) == level(s) #=
    =# && spin(p) != spin(q) && spin(r) != spin(s)
        spin(p) == spin(r) ? -g/2 : g/2
    else
        0
    end
end

function V_mb(g)
#    V = V_sp(g)
    @Operator(MBBASIS) do X, Y
        println("V: ", (index(X), index(Y)))

        sum(cartesian_pow(SPBASIS, Val{2})) do I
            p, r = I
            q = flipspin(p)
            s = flipspin(r)

            -g/2*A(p', q', s, r)(X, Y)
        end
    end
end

H(g) = @Operator(MBBASIS) do X, Y
    println("H: ", (index(X), index(Y)))
    f_mb(X, Y) + V_mb(g)(X, Y)
end

main() = tabulate(H(1.0))

end # module Exec

Exec.main()
