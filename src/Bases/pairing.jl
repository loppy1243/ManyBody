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
flipspin(p::Sub{<:Pairing}) = typeof(p)(flipspin(p))
flipspin(p::Index{<:Pairing}) = typeof(p)(flipspin(p))

Base.:(==)(p1::Pairing{L}, p2::Pairing{L}) where L =
    p1.level == p2.level && p1.spin == p2.spin

index(sp::Pairing) = 2(sp.level-1) + 1 + Bool(sp.spin)

index_to_spin(i) = Spin(Bool(1 - i % 2))
index_to_level(i) = div(i - Bool(index_to_spin(i)) - 1, 2) + 1

indexbasis(::Type{Pairing{L}}, pn::Int) where L =
    Pairing{L}(index_to_spin(pn), index_to_level(i))

nlevels(::Type{Pairing{L}}) where L = L
nlevels(::Type{<:Sub{Pairing{L}}}) where L = L
nlevels(::Type{Index{Pairing{L}}}) where L = L
level(sp::Pairing) = sp.level
level(sp::Sub{<:Pairing}) = sp.state.level
level(sp::Index{<:Pairing}) = index_to_level(sp.index)
spin(sp::Pairing) = sp.spin
spin(sp::Sub{<:Pairing}) = sp.state.spin
spin(sp::Index{<:Pairing}) = index_to_level(sp.index)

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
