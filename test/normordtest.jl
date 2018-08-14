function normordtest()
@testset "Normal Ordering" begin
    SPBASIS = Bases.Pairing{4}
    p = SPBASIS(1, SPINUP)
    q = SPBASIS(2, SPINDOWN)
    r = SPBASIS(4, SPINUP)
    s = SPBASIS(4, SPINDOWN)
    a = A(p', r, s', q')

    @testset "Vacuum" begin
        @debug testset="Normal Ordering/Vacuum"
        sgn, Na = normord(a)
        @debug begin
            input   = "+ $a"
            want    = "+ A($s', $q', $p', $r )"
            normord = "$(sgn > 0 ? "+" : "-") $Na"
        end
        @test Na == A(s', q', p', r)
        @test sgn == 1
    end

    @testset "Fermi{2}" begin
        @debug testset="Normal Ordering/Fermi{2}"
        sgn, Na = normord(RefStates.Fermi{2, SPBASIS}, a)
        @debug begin
            input   = "+ $a"
            want    = "- A($s', $p', $q', $r )"
            normord = "$(sgn > 0 ? "+" : "-") $Na"
        end
        @test Na == A(s', p', q', r)
        @test sgn == -1
    end
end; nothing end
normordtest()
