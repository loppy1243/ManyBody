module Exec

include("src/ManyBody.jl")

using Loppy.Util: cartesian_pow
using .ManyBody

const LEVEL_SPACING = 1
const SPBASIS = Bases.Pairing{4}
const REFSTATE = RefStates.Fermi{2, SPBASIS}
const MBBASIS = Bases.Paired{2, 4}
#const INDICES = indices(SPBASIS)
#const PARTS = parts(REFSTATE)
#const HOLES = holes(REFSTATE)

E0(g) = sum(SPBASIS) do p
    isocc(REFSTATE, p)*(LEVEL_SPACING*(level(p)-1) - g/2)
end

f(g) = @Operator(MBBASIS) do X, Y
    sum(SPBASIS) do p
#        sgn1, NA = normord(REFSTATE, A(p', p))
        sgn1, NA = 1, A(p', p)
        sgn2, Y2 = apply_normord_rl(NA, Y)
        (LEVEL_SPACING*(level(p)-1))*sgn1*sgn2*(X'Y2)
    end
end

V(g) = @Operator(MBBASIS) do X, Y
    X == Y ? -g : index(X) + index(Y) - 1 == dim(MBBASIS) ? 0.0 : -g/2
#    sum(cartesian_pow(SPBASIS, Val{2})) do I
#        p, r = I
#        q = flipspin(p)
#        s = flipspin(r)
##        sgn1, NA = normord(REFSTATE, a)
#        sgn1, NA = 1, A(p', q', s, r)
#        sgn2, Y2 = apply_normord_rl(NA, Y)
#
#        -g/2*sgn1*sgn2*(X'Y2)
#    end
end

H(g) = @Operator(MBBASIS) do X, Y
    #=(X'Y)*E0(g) +=# f(g)(X, Y) + V(g)(X, Y)
end

main() = tabulate(H(1.0)).op

end # module Exec
