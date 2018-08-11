computed_default(BASIS) = @testset "Computed Default" begin
    @test Bases.Generation(BASIS) === Bases.Computed()
    @test basis(BASIS) !== basis(BASIS)
end

generate_macro(BASIS) = @testset "Bases.@generate" begin
    @test Bases.Generation(BASIS) === Bases.Generated()
    @test basis(BASIS) === basis(BASIS)
end

@testset "Basis Generation" begin
    BASIS = BASES.Pairing{4}

    computed_default(BASIS)
    Bases.@generate BASIS
    generate_macro(BASIS)
end
