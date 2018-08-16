productbasis() = @testset "Product Basis" begin
    B = Product{3, NTuple{3, Bases.Pairing{4}}}
    b = B(Bases.Pairing{4}.((1, 1, 2, 2), (SPINDOWN, SPINUP, SPINDOWN, SPINUP)))

    i = index(b)
    b2 = B[index(b)]
    @debug "Product Basis" index=i element=b2

    @test b == b2
end
productbasis()
