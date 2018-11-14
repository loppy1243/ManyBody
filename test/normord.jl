using ManyBody
using ManyBody.Operators: @A

normordtest() = @testset "Normal Ordering" begin
    SPBASIS = Bases.Pairing{4}
    MBBASIS = Bases.Slater{SPBASIS}
    p = SPBASIS(1, SPINUP)
    q = SPBASIS(2, SPINDOWN)
    r = SPBASIS(4, SPINUP)
    s = SPBASIS(4, SPINDOWN)
    a = @A(p', r, s', q')

    @testset "Vacuum" begin
        @debug "Testing normord wrt. Vacuum"

        Na_mat = tabulate(normord(a), Int, MBBASIS)
        correct_mat = tabulate(@A(s', q', p', r), Int, MBBASIS)
        @test all(Na_mat .== correct_mat)
    end

    @testset "Fermi{2}" begin
        @debug "Testing normord wrt. Fermi{2}",

        Na_mat = tabulate(normord(RefStates.Fermi{SPBASIS}(2), a), Int, MBBASIS)
        correct_mat = tabulate(@A(s', p', q', r), Int, MBBasis)
        @test all(Na_mat .== correct_mat)
    end
end
