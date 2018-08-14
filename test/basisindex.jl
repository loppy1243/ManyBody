function basisindex()
@testset "Basis Indexing" begin
    for BASIS in (Bases.Pairing{4},
                  Bases.Slater{Bases.Pairing{4}},
                  Bases.Paired{2, 4})
        b = basis(BASIS)

        @debug "Testing basis: " BASIS
        @test length(b) == dim(BASIS)
        @test all(b[i] == BASIS[i] for i in 1:dim(BASIS))
    end
end; end
basisindex()
