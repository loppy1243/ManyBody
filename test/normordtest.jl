function normordtest()
@testset "Normal Ordering" begin
    SPBASIS = Bases.Pairing{4}
    p = SPBASIS(1, SPINUP)
    q = SPBASIS(2, SPINDOWN)
    r = SPBASIS(4, SPINUP)
    s = SPBASIS(4, SPINDOWN)
    a = A(p', r, s', q')

    @testset "Vacuum" begin
        sgn, Na = normord(a)
        @debug("Testing normord wrt. Vacuum",
               input   = (1, a),
               want    = (1, A(s', q', p', r)),
               normord = (sgn, Na))
        @test Na == A(s', q', p', r)
        @test sgn == 1
    end

    @testset "Fermi{2}" begin
        sgn, Na = normord(RefStates.Fermi(SPBASIS, 2), a)
        @debug("Testing normord wrt. Fermi{2}",
               input   = (1, a),
               want    = (-1, A(s', p', q', r)),
               normord = (sgn, Na))
        @test Na == A(s', p', q', r)
        @test sgn == -1
    end
end; nothing end
normordtest()
