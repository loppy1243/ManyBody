export Pairing, level, nlevels, spin, flipspin
import .RefStates
using ..SpinMod

struct Pairing{Levels} <: AbstractBasis
    level::Int
    spin::Spin

    function Pairing{Levels}(l::Int, s::Spin) where Levels
        @assert l <= Levels
        new(l, s)
    end
end

flipspin(p::Pairing) = typeof(p)(level(p), flip(spin(p)))

Base.:(==)(p1::Pairing{L}, p2::Pairing{L}) where L =
    p1.level == p2.level && p1.spin == p2.spin

index(sp::Pairing) = 2(sp.level-1) + 1 + Bool(sp.spin)

function indexbasis(::Type{Pairing{L}}, pn::Int) where L
    spin_num = Bool(1 - pn % 2)
    level = div(pn - spin_num - 1, 2) + 1

    Pairing{L}(level, Spin(spin_num))
end

nlevels(::Type{Pairing{L}}) where L = L
level(sp::Pairing) = sp.level
spin(sp::Pairing) = sp.spin

dim(::Type{Pairing{L}}) where L = 2L

Base.show(io::IO, x::Pairing) = print(io, level(x), spinup(spin(x)) ? "↑" : "↓")
function Base.show(io::IO, ::MIME"text/plain", x::Pairing{L}) where L
    for i = L:-1:1
        print(io, lpad(i, ndigits(L)), "|  ")
        if i == level(x)
            if spindown(spin(x))
                print(io, "↓  ")
            else
                print(io, "o  ")
            end

            if spinup(spin(x))
                println(io, "↑")
            else
                println(io, "o")
            end
        else
            println(io, "o  o")
        end
    end
end
