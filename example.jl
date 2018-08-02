module Exec

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

E0(g) = sum(SPBASIS) do p
    isocc(REFSTATE, p)*(LEVEL_SPACING*(level(p)-1) - g/2)
end

f(g) = @Operator(MBBASIS) do X, Y
    sum(SPBASIS) do p
        sgn1, NA = normord(REFSTATE, A(p', p))
        sgn2, Y2 = apply_normord_rl(NA, Y)
        (LEVEL_SPACING*(level(p)-1) - g*isocc(REFSTATE, p))*sgn1*sgn2*(X'Y2)
    end
end

#V_sp(g) = @Operator(SPBASIS) do p, q, r, s
#    if level(p) == level(q) && level(r) == level(s) #=
#    =# && spin(p) != spin(q) && spin(r) != spin(s)
#        spin(p) == spin(r) ? -g/2 : g/2
#    else
#        0
#    end
#end

V(g) = @Operator(MBBASIS) do X, Y
    sum(cartesian_pow(SPBASIS, Val{2})) do I
        p, r = I
        q = flipspin(p)
        s = flipspin(r)
        sgn1, NA = normord(REFSTATE, A(p', q', s, r))
        sgn2, Y2 = apply_normord_rl(NA, Y)

        -g/2*sgn1*sgn2*(X'Y2)
    end
end

H(g) = @Operator(MBBASIS) do X, Y
    (X'Y)*E0(g) + f(g)(X, Y) + V(g)(X, Y)
end

main() = tabulate(H(1.0))

end # module Exec

#Exec.main()
