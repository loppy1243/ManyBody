using ManyBody

const LEVELS = 4
const BASIS = States.PairingParticle{LEVELS}
const REFSTATE = RefStates.Fermi{2, BASIS}
const PARTS = particles(REFSTATE)
const HOLES = holes(REFSTATE)

V(g) = Operator{2, BASIS, Float64}() do p, q, r, s
    if level(p) == level(q) && level(r) == level(s) #=
    =# && spin(p) != spin(q) && spin(r) != spin(s)
        spin(p) == spin(r) ? -g/2 : g/2
    else
        0
    end
end |> tabulate
