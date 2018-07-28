export PairingParticle, level, snum, particle, nlevels, spin, states, nstates

import .RefStates
using ..SpinMod

struct PairingParticle{Levels} <: SPState
    level::Int
    spin::Spin

    function PairingParticle{Levels}(l::Int, s::Spin) where Levels
        @assert l <= Levels
        new(l, s)
    end
end

Base.:(==)(::PairingParticle, ::PairingParticle) = false
Base.:(==)(p1::PairingParticle{L}, p2::PairingParticle{L}) where L =
    p1.level == p2.level && p1.spin == p2.spin

snum(sp::Int) = sp
snum(sp::PairingParticle) = 2(s.level-1) + 1 + Bool(s.spin)

function particle(::Type{PairingParticle{L}}, pn::Int) where L
    spin_num = Bool(1 - pn % 2)
    level = div(pn - spin - 1, 2) + 1

    PairingParticle{L}(level, Spin(spin_num))
end

nlevels(::Type{PairingParticle{L}}) where L = L
level(sp::PairingParticle) = s.level
spin(sp::PairingParticle) = s.spin

states(::Type{PairingParticle{L}}) where L =
    (PairingParticle{L}(l, s) for l = 1:L, s in SPINS)

nstates(p::PairingParticle{L}) where L = 2L
