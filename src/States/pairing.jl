struct PairingParticle{Levels} <: SPState
    level::Int
    spin::Spin

    function PairingParticle(l::Int, s::Spin)
        @assert l <= Levels
        new(l, s)
    end
end

pnum(sp::Int) = sp
pnum(sp::PairingParticle) = 2(s.level-1) + 1 + Bool(s.spin)

function particle(::Type{PairingParticle{L}}, pn::Int) where L
    spin_num = Bool(1 - pn % 2)
    level = div(pn - spin - 1, 2) + 1

    PairingParticle{L}(level, Spin(spin_num))
end

nlevels(sp::PairingParticle{L}) where L = L
level(sp::PairingParticle) = s.level
spin(sp::PairingParticle) = s.spin

iter(P::Type{PairingParticle{L}}) where L = (PairingParticle{L}(l, s) for l = 1:L, s in SPINS)
