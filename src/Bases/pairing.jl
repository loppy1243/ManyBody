export Pairing, level, nlevels, spin, flipspin
import .RefStates
using ..SpinMod

struct Pairing{Levels} <: ConcreteBasis
    level::Int
    spin::Spin

    function Pairing{Levels}(l::Int, s::Spin) where Levels
        @assert l <= Levels
        new(l, s)
    end
end
IndexType(::Type{<:Pairing}) = IndexTypes.Linear()

flipspin(p::Pairing) = typeof(p)(level(p), flip(spin(p)))
flipspin(p::Wrapped) = convert(typeof(p), flipspin(inner(p)))

Base.:(==)(p1::Pairing{L}, p2::Pairing{L}) where L =
    p1.level == p2.level && p1.spin == p2.spin

index(sp::Pairing) = 2(sp.level-1) + 1 + Bool(sp.spin)

function indexbasis(B::Type{<:Pairing}, pn::Int) where L
    s_num = Bool(1 - i % 2)
    l = div(i - s_num - 1, 2) + 1
    B(l, Spin(s_num))
end

nlevels(a::AbstractBasis) = nlevels(typeof(a))
nlevels(::Type{Pairing{L}}) where L = L
nlevels(B::Type{<:Wrapped}) = nlevels(innertype(B))

level(sp::Pairing) = sp.level
level(sp::Wrapped) = level(inner(sp))

spin(sp::Pairing) = sp.spin
spin(sp::Wrapped) = spin(inner(sp))

SpinMod.spinup(sp) = spinup(spin(sp))

dim(B::Type{<:Pairing}) = 2nlevel(B)

Base.show(io::IO, x::Pairing) = print(io, level(x), spinup(spin(x)) ? "↑" : "↓")
function Base.show(io::IO, ::MIME"text/plain", x::Pairing)
    L = nlevels(x)
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
