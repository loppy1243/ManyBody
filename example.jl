using ManyBody

const LEVELS = 4
const BASIS = States.PairingParticle{LEVELS}
const REFSTATE = RefStates.Fermi{2, BASIS}
const PARTS = particles(REFSTATE)
const HOLES = holes(REFSTATE)

const f_spmat = Operator1Body{BASIS}() do p, q
    LEVEL_SPACING*(level(p)-1)
end |> tabulate

const f_mb = MBOperator{MBState{REFSTATE, IntState{BASIS}}}() do X, Y
    ssi = intstates(BASIS)

    sum(ssi) do p
        f_spmat(p, p)*A(p', p)(X, T)
    end
end

V_sp(g) = Operator2Body{BASIS}() do p, q, r, s
    if level(p) == level(q) && level(r) == level(s) #=
    =# && spin(p) != spin(q) && spin(r) != spin(s)
        spin(p) == spin(r) ? -g/2 : g/2
    else
        0
    end
end |> tabulate

function V_mb(g)
    V = V_sp(g)
    MBOperator{MBState{REFSTATE, IntState{BASIS}}}() do X, Y
        sum(Iterators.product(fill(intstates(BASIS), 4))) do (p, q, r, s)
            V(p, q, r, s)*A(p', q', s, r)(X, Y)
        end
    end
end

H(g) = MBOperator{MBState{REFSTATE, IntState{BASIS}}}() do X, Y
#    zerobody = ...
    zerobody + f_mb(X, Y) + V_mb(g)(X, Y)
end
