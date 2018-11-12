module Exec

using ManyBody
using ManyBody.Operators: @A

Bases.@defSub(NPairing{L, P} <: Bases.Slater{Bases.Pairing{L}}) do b
    nocc(b) == P && all(p -> p.level <= L, occ(b))
end

const LEVEL_SPACING = 1
const SPBASIS = Bases.Pairing{4}
const REFSTATE = RefStates.Fermi{SPBASIS, 2}
const MBBASIS = Bases.Paired{2, 4}
#const MBBASIS = NPairing{4, 4}

f(g) = (X, Y) -> sum(spbasis(Y)) do p
    sgn, Y′ = @A(p', p)
    LEVEL_SPACING*(p.level-1)*sgn*overlap(X, Y′)
end

V(g) = (X, Y) -> begin
    SPB = spbasis(Y)
    -g/2 * sum(Iterators.product(SPB, SPB, SPB, SPB)) do (p, q, r, s)
        mask  = p.level == q.level
        mask *= r.level == s.level
        mask *= spinup(p)*spindown(q)
        mask *= spinup(r)*spindown(s)

        sgn, Y′ = @A(p', q', s, r)
        mask*sgn*overlap(X, Y′)
    end
end

H(g) = f(g) + V(g)
main() = tabulate(H(1.0), Float64, MBBASIS, MBBASIS)

end # module Exec
Exec.main()
