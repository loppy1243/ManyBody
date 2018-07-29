module Exec

include("src/ManyBody.jl")

using Loppy.Util: cartesian_pow
using .ManyBody

const LEVEL_SPACING = 1
@default_basis! BASIS = Bases.Pairing{4}
@default_refstate! REFSTATE = RefStates.Fermi{2, BASIS}
const MBBASIS = MBBasis{REFSTATE, BASIS}
#const INDICES = indices(BASIS)
#const PARTS = parts(REFSTATE)
#const HOLES = holes(REFSTATE)

const f_sp = @Operator() do p, q
    LEVEL_SPACING*(level(p)-1)
end

const f_mb = @MBOperator() do X, Y
    println("f: ", (index(X), index(Y)))
    sum(BASIS) do p
        f_sp(p, p)*A(p', p)(X, Y)
    end
end

V_sp(g) = @Operator() do p, q, r, s
    if level(p) == level(q) && level(r) == level(s) #=
    =# && spin(p) != spin(q) && spin(r) != spin(s)
        spin(p) == spin(r) ? -g/2 : g/2
    else
        0
    end
end

function V_mb(g)
    V = V_sp(g)
    @MBOperator() do X, Y
        println("V: ", (index(X), index(Y)))
        sum(cartesian_pow(BASIS, Val{2})) do I
            p, r = I
            q = flipspin(p)
            s = flipspin(r)

            V(p, q, r, s)*A(p', q', s, r)(X, Y)
        end
    end
end

H(g) = @MBOperator() do X, Y
    println("H: ", (index(X), index(Y)))
    f_mb(X, Y) + V_mb(g)(X, Y)
end

main() = tabulate(H(1.0))

end # module Exec

Exec.main()
