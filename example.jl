using ManyBody

const LEVELS = 4
const BASIS = States.PairingParticle{LEVELS}
const REFSTATE = RefStates.Fermi{2, BASIS}
const PARTS = particles(REFSTATE)
const HOLES = holes(REFSTATE)

V_sp(g) = SPOperator{2, BASIS, Float64}() do p, q, r, s
    if level(p) == level(q) && level(r) == level(s) #=
    =# && spin(p) != spin(q) && spin(r) != spin(s)
        spin(p) == spin(r) ? -g/2 : g/2
    else
        0
    end
end |> tabulate

V_mb(g) = MBOperator{MBState{REFSTATE, BASIS}, Float64}() do X, Y
    sum(Iterators.product(fill(states(BASIS), 4))) do (p, q, r, s)
        V_sp(g)*A(p', q', s, r)(X, Y)
    end
end

H(g) = MBOperator{MBState{REFSTATE, BASIS}, Float64}() do X, Y
    ss = states(BASIS)

    zerobody = 0

    onebody = sum(ss) do p
        LEVEL_SPACING*(level(p)-1)*A(p', p)(X, Y)
    end

    twobody = sum(Iterators.product(ss, ss, ss, ss)) do (p, q, r, s)
        V_sp(g)(p, q, r, s)*A(p', q', s, r)(X, Y)
    end
end
