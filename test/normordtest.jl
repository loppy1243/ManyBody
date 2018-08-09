function normordtest()
@testset "Normal Ordering" begin
    SPBASIS = Bases.Pairing{4}
    p = SPBASIS(1, SPINUP)
    q = SPBASIS(2, SPINDOWN)
    r = SPBASIS(4, SPINUP)
    s = SPBASIS(4, SPINDOWN)
    a = A(p', r, s', q')

    @testset "Vacuum" begin
        println("Vacuum")
        sgn, Na = normord(a)
        println("Input:")
        println("+ $a")
        println("Want:")
        println("+ A($s', $q', $p', $r )")
        println("Normal ordered:")
        print(sgn > 0 ? "+ " : "- ", Na)
        println()
        println()
        @test Na == A(s', q', p', r)
        @test sgn == 1
    end

    @testset "Fermi{2}" begin
        println("Fermi{2}")
        sgn, Na = normord(RefStates.Fermi{2, SPBASIS}, a)
        println("Input:")
        println("+ $a")
        println("Want:")
        println("- A($s', $p', $q', $r )")
        println("Normal ordered:")
        println(sgn > 0 ? "+ " : "- ", Na)
        println()
        println()
        @test Na == A(s', p', q', r)
        @test sgn == -1
    end
end; nothing end
normordtest()
