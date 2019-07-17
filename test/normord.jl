using ManyBody
using ManyBody.Operators: @A

normordtest() = @testset "Normal Ordering" begin
    SPBASIS = Bases.Pairing{4}
    MBBASIS = Bases.Slater{SPBASIS, rank(SPBASIS)}
    p = SPBASIS(1, SPINUP)
    q = SPBASIS(2, SPINDOWN)
    r = SPBASIS(4, SPINUP)
    s = SPBASIS(4, SPINDOWN)
    a = @A(p', r, s', q')

    @testset "Vacuum" begin
        @debug "Testing normord wrt. Vacuum"

        sgn, Na = normord(a)
        Na_mat = sgn.*Na.(MBBASIS)
        correct_mat = @A(s', q', p', r).(MBBASIS)
        @test all(Na_mat .== correct_mat)
    end

    @testset "Fermi{2}" begin
        @debug "Testing normord wrt. Fermi{2}"

        sgn, Na = normord(RefStates.Fermi{SPBASIS}(2), a)
        Na_mat = sgn.*Na.(MBBASIS)
        correct_mat = .-@A(s', r, p', q').(MBBASIS)
        @test all(Na_mat .== correct_mat)
    end
end
