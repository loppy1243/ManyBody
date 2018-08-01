using BenchmarkTools: @btime

function rlexpect()
    SPBASIS = Bases.Pairing{4}
    REFSTATE = RefStates.Fermi{2, SPBASIS}
    MBBASIS = Bases.PartHole{RefStates.Fermi{2, SPBASIS}}
    p, q, r, s = rand(basis(SPBASIS), 4)
    h, i = shuffle(parts(REFSTATE))[1:2]
    j, k = shuffle(parts(REFSTATE))[1:2]
    l, m = shuffle(holes(REFSTATE))[1:2]
    n, o = shuffle(holes(REFSTATE))[1:2]
    X = MBBASIS((h, l), (i, m))
    Y = MBBASIS((j, n), (k, o))
    a = A(p', q', r, s)

    println("Contraction method:")
    @btime $a($X, $Y)
    println("Direct method:")
    @btime ($X')*($a($Y))

    nothing
end
