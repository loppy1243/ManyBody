rlapplytest() = @testset "RaiseLowerOp Application" begin
    SPBASIS = Bases.Pairing{4}
    MBBASIS = Bases.Paired{2, 4}

    p, q, r, s, t, u = SPBASIS.((1, 3, 2, 3, 4, 1), (↑, ↓, ↑, ↑, ↓, ↓))

    ==′((sgn1, X′), (sgn2, X)) = sgn1 == sgn2 && (iszero(sgn1) || X′ == X)
    @testset "One Body" begin for X in +MBBASIS
        X = convert(supbasis(MBBASIS), X)
        @test @A(p, p)(X)[1] |> iszero
        @test @A(q, q)(X)[1] |> iszero
        @test @A(p', p')(X)[1] |> iszero
        @test @A(q', q')(X)[1] |> iszero
        @test @A(p', p)(X) ==′ (p in X, X)
        @test @A(q', q)(X) ==′ (q in X, X)
        @test @A(p, p')(X) ==′ (1 - (p in X), X)
        @test @A(q, q')(X) ==′ (1 - (q in X), X)

        if p != q && !(p in X) || !(q in X)
            @test @A(p, q)(X)[1] |> iszero
            @test @A(q, p)(X)[1] |> iszero
        end
        if p != q && p in X || q in X
            @test @A(p', q')(X)[1] |> iszero
            @test @A(q', p')(X)[1] |> iszero
        end
    end end

#    ## Need to rewrite these.
#    @testset "Two Body" begin
#        Z = convert(supbasis(MBBASIS), (-MBBASIS)[1])
#        @test @A(p, q, r, s)(Z)[1] |> iszero
#        X = Bases.Slater(SPBASIS.((1, 2, 3, 3),
#                                  (SPINDOWN, SPINDOWN, SPINDOWN, SPINUP)))
#        sgn, Y = @A(p, q', r, s')(Z)
#        @test sgn == -1 && X == Y
#        X = Bases.Slater(SPBASIS.((1, 2), (SPINDOWN, SPINDOWN)))
#        sgn, Y = @A(p, q, r, q')(Z)
#        @test sgn == -1 && X == Y
#        X = Bases.Slater(SPBASIS.((1, 2, 3, 4), (SPINUP, SPINDOWN, SPINDOWN, SPINDOWN)))
#        sgn, Y = @A(t', r, u, q')(Z) 
#        @test sgn == -1 && X == Y
#    end
end
