using BenchmarkTools: @btime

function rlexpect()
    SPBASIS = Bases.Pairing{4}
    REFSTATE = RefStates.Fermi{2, SPBASIS}
    PARTS = parts(REFSTATE)
    HOLES = holes(REFSTATE)
    MBBASIS = Bases.Slater{SPBASIS}

    p, q, r, s = rand(basis(SPBASIS), 4)
    h, i = shuffle(PARTS)[1:2]
    j, k = shuffle(PARTS)[1:2]
    l, m = shuffle(HOLES)[1:2]
    n, o = shuffle(HOLES)[1:2]

    X = MBBASIS(h, l, i, m)
    Y = MBBASIS(j, n, k, o)

    a = A(p', q', r, s)

    println("Direct method:")
    @btime ($X')*($a($Y))

    nothing
end
