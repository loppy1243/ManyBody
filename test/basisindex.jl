function basisindex()
@testset "Basis Indexing" begin
    for BASIS in (Bases.Pairing{4},
                  Bases.PartHole{RefStates.Fermi{2, Bases.Pairing{4}}},
                  Bases.PHPaired{2, 4})
        b = basis(BASIS)

        println("Testing basis: ", BASIS)
        @test length(b) == dim(BASIS)
        @test all(b[i] == BASIS[i] for i in indices(BASIS))
    end
end; end
