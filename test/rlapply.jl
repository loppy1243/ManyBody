function rlapply()
@testset "RaiseLowerOp Application" begin
    SPBASIS = Bases.Pairing{4}
#    REFSTATE = RefStates.Fermi{2}
    MBBASIS = Bases.Paired{2, 4}

    p, q = SPBASIS.((1, 3), (SPINUP, SPINDOWN))

    for X in MBBASIS
        @test apply_normord_rl(A(p, p), X) == (0, States.ZERO)
        @test apply_normord_rl(A(q, q), X) == (0, States.ZERO)
        @test apply_normord_rl(A(p', p'), X) == (0, States.ZERO)
        @test apply_normord_rl(A(q', q'), X) == (0, States.ZERO)
        @test apply_normord_rl(A(p', p), X) == (p in X, p in X ? X : States.ZERO)
        @test apply_normord_rl(A(q', q), X) == (q in X, q in X ? X : States.ZERO)
        @test apply_normord_rl(A(p, p'), X) == (-(p in X), p in X ? X : States.ZERO)
        @test apply_normord_rl(A(q, q'), X) == (-(q in X), q in X ? X : States.ZERO)

        if !(p in X) || !(q in X)
            @test apply_normord_rl(A(p, q), X) == (0, States.ZERO)
            @test apply_normord_rl(A(q, p), X) == (0, States.ZERO)
        end
        if p in X || q in X
            @test apply_normord_rl(A(p', q'), X) == (0, States.ZERO)
            @test apply_normord_rl(A(q', p'), X) == (0, States.ZERO)
        end
    end
end; end
