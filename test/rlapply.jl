function rlapply()
@testset "RaiseLowerOp Application" begin
    SPBASIS = Bases.Pairing{4}
    MBBASIS = Bases.Paired{2, 4}

    p, q, r, s = SPBASIS.((1, 3, 2, 3), (SPINUP, SPINDOWN, SPINUP, SPINUP))

    @testset "One Body" begin for X in MBBASIS
        X2 = convert(Bases.Slater, X)
        @test A(p, p)(X) == (0, States.ZERO)
        @test A(q, q)(X) == (0, States.ZERO)
        @test A(p', p')(X) == (0, States.ZERO)
        @test A(q', q')(X) == (0, States.ZERO)
        @test A(p', p)(X) == (p in X, p in X ? X2 : States.ZERO)
        @test A(q', q)(X) == (q in X, q in X ? X2 : States.ZERO)
        @test A(p, p')(X) == (1 - (p in X), p in X ? States.ZERO : X2)
        @test A(q, q')(X) == (1 - (q in X), q in X ? States.ZERO : X2)

        if p != q && !(p in X) || !(q in X)
            @test A(p, q)(X) == (0, States.ZERO)
            @test A(q, p)(X) == (0, States.ZERO)
        end
        if p != q && p in X || q in X
            @test A(p', q')(X) == (0, States.ZERO)
            @test A(q', p')(X) == (0, States.ZERO)
        end
    end end

    @testset "Two Body" begin
        X = Bases.Slater(SPBASIS.((1, 2, 3, 3),
                                  (SPINDOWN, SPINDOWN, SPINDOWN, SPINUP)))
        @test A(p, q, r, s)(MBBASIS[1]) == (0, States.ZERO)
        @test A(p, q', r, s')(MBBASIS[1]) #=
           =# == (-1, X)
    end
end; nothing end
