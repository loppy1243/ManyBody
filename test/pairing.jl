using ManyBody
using ManyBody.Operators: @A

#Bases.@defSub NPairing{L, P} <: Bases.Slater{Bases.Pairing{L}} begin b
#    nocc(b) == P && all(p -> p.level <= L, occ(b))
#end

_pairingH_direct(g, i, j) = if i == j
    [2, 4, 6, 6, 8, 10][i]*LEVEL_SPACING - g
elseif i + j - 1 == 6
    0.0
else
    -G/2
end
pairingH_direct(g) = [_H(g, i, j) for i = 1:6, j = 1:6]

pairingtest(; g_samples, atol) = @testset "Pairing Hamiltonian" begin
    const LEVEL_SPACING = 1
    const SPBASIS = Bases.Pairing{4}
    const REFSTATE = RefStates.Fermi{SPBASIS}(2)
    const MBBASIS = Bases.Paired{2, 4}
    
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
    
    H(g) = (X, Y) -> f(g)(X, Y) + V(g)(X, Y)

    for g in range(-1.0, stop=1.0, length=5)
        @test all(abs.(tabulate(H(g), Float64, MBBASIS) - pairingH_direct(g)) #=
               =# .< atol)
    end
end
