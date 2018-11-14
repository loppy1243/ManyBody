rlapplytest() = @testset "RaiseLowerOp Application" begin
    SPBASIS = Bases.Pairing{4}
    MBBASIS = Bases.Paired{2, 4}

    p, q, r, s, t, u = SPBASIS.((1, 3, 2, 3, 4, 1),
                                (SPINUP, SPINDOWN, SPINUP, SPINUP, SPINDOWN, SPINDOWN))

    @testset "One Body" begin for X in MBBASIS
        @test @A(p, p)(X)[1] == 0
        @test @A(q, q)(X)[1] == 0
        @test @A(p', p')(X)[1] == 0
        @test @A(q', q')(X)[1] == 0
        @test @A(p', p)(X) == (p in X, X)
        @test @A(q', q)(X) == (q in X, X)
        @test @A(p, p')(X) == (1 - (p in X), X)
        @test @A(q, q')(X) == (1 - (q in X), X)

        if p != q && !(p in X) || !(q in X)
            @test A(p, q)(X)[1] == 0
            @test A(q, p)(X)[1] == 0
        end
        if p != q && p in X || q in X
            @test A(p', q')(X)[1] == 0
            @test A(q', p')(X)[1] == 0
        end
    end end

    @testset "Two Body" begin
        @test @A(p, q, r, s)(MBBASIS[1])[1] == 0
        X = Bases.Slater(SPBASIS.((1, 2, 3, 3),
                                  (SPINDOWN, SPINDOWN, SPINDOWN, SPINUP)))
        sgn, Y = @A(p, q', r, s')(MBBASIS[1])
        @test sgn == -1 && X == Y
        X = Bases.Slater(SPBASIS.((1, 2), (SPINDOWN, SPINDOWN)))
        sgn, Y = @A(p, q, r, q')(MBBASIS[1])
        @test sgn == -1 && X == Y
        X = Bases.Slater(SPBASIS.((1, 2, 3, 4), (SPINUP, SPINDOWN, SPINDOWN, SPINDOWN)))
        sgn, Y = A(t', r, u, q')(MBBASIS[1]) 
        @test sgn == -1 && X == Y
    end
end
