function normordtest()
@testset "Normal Ordering" begin
    SPBASIS = Bases.Pairing{4}
    p = SPBASIS(1, SPINUP)
    q = SPBASIS(2, SPINDOWN)
    r = SPBASIS(4, SPINUP)
    s = SPBASIS(4, SPINDOWN)
    a = A(p', r, s', q')

    sgn, Na = normord(a)
    println("Input:")
    map(x -> print(x.state, " "), a.ops)
    println()
    println("Want:")
    println("$s $q $p $r")
    println("Normal ordered:")
    map(x -> print(x.state, " "), Na.ops)
    println()
    println()
    @test Na == A(s', q', p', r)
    @test sgn == 1
end; end
