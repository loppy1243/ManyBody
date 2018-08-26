function rlapply()
@testset "RaiseLowerOp Application" begin
    SPBASIS = Bases.Pairing{4}
    MBBASIS = Bases.Paired{2, 4}

    p, q, r, s, t, u = SPBASIS.((1, 3, 2, 3, 4, 1),
                                (SPINUP, SPINDOWN, SPINUP, SPINUP, SPINDOWN, SPINDOWN))

    @testset "One Body" begin for X in MBBASIS
        @test A(p, p)(X) == States.Zero()
        @test A(q, q)(X) == States.Zero()
        @test A(p', p')(X) == States.Zero()
        @test A(q', q')(X) == States.Zero()
        @test A(p', p)(X) == (p in X) * X
        @test A(q', q)(X) == (q in X) * X
        @test A(p, p')(X) == (1 - (p in X)) * X
        @test A(q, q')(X) == (1 - (q in X)) * X

        if p != q && !(p in X) || !(q in X)
            @test A(p, q)(X) == States.Zero()
            @test A(q, p)(X) == States.Zero()
        end
        if p != q && p in X || q in X
            @test A(p', q')(X) == States.Zero()
            @test A(q', p')(X) == States.Zero()
        end
    end end

    @testset "Two Body" begin
        @test A(p, q, r, s)(MBBASIS[1]) == States.Zero()
        X = Bases.Slater(SPBASIS.((1, 2, 3, 3),
                                  (SPINDOWN, SPINDOWN, SPINDOWN, SPINUP)))
        @test A(p, q', r, s')(MBBASIS[1]) == -X
        X = Bases.Slater(SPBASIS.((1, 2), (SPINDOWN, SPINDOWN)))
        @test A(p, q, r, q')(MBBASIS[1]) == -X
        X = Bases.Slater(SPBASIS.((1, 2, 3, 4), (SPINUP, SPINDOWN, SPINDOWN, SPINDOWN)))
        @test A(t', r, u, q')(MBBASIS[1]) == -X
    end
end; nothing end
rlapply()
