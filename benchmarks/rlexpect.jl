using BenchmarkTools: @btime

function rlexpect()
    SPBASIS = Bases.Pairing{4}
    MBBASIS = Bases.PartHole{RefStates.Fermi{2, SPBASIS}}
    X, Y = rand(basis(MBBASIS), 2)
    p, q, r, s = rand(basis(SPBASIS), 4)

    println("Contraction method:")
    @btime $(A(p', q', r, s))(X, Y)
    println("Direct method:")
    @btime X'$(A(p', q', r, s))(Y)

    nothing
end
