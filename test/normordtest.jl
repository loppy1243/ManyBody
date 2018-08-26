function normordtest()
@testset "Normal Ordering" begin
    SPBASIS = Bases.Pairing{4}
    p = SPBASIS(1, SPINUP)
    q = SPBASIS(2, SPINDOWN)
    r = SPBASIS(4, SPINUP)
    s = SPBASIS(4, SPINDOWN)
    a = A(p', r, s', q')

    @testset "Vacuum" begin
        @debug "Testing normord wrt. Vacuum"

        Na = normord(a)
        @test tabulate(Na).rep == tabulate(A(s', q', p', r)).rep
    end

    @testset "Fermi{2}" begin
        @debug "Testing normord wrt. Fermi{2}",

        Na = normord(RefStates.Fermi(SPBASIS, 2), a)
        @test tabulate(Na).rep == tabulate(A(s', p', q', r)).rep
    end
end; nothing end
normordtest()
