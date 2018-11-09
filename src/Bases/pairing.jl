export Pairing, nlevels, flipspin
#import .RefStates
using ..SpinMod

SpinMod.spinup(sp::Pairing) = spinup(sp.spin)

struct Pairing{Levels} <: TensorBasis{2}
    level::Int
    spin::Spin

    function Pairing{L}(l::Int, s::Spin) where L
        @assert l <= L
        new(l, s)
    end
end

fulldims(B::Type{<:Pairing}) = (nlevels(B), 2)
index(sp::Pairing) = CartesianIndex(p.level, 1+Bool(p.spin))
indexbasis(B::Type{<:Pairing}, level::Int, snum::Int) =
    B(level, Spin(Bool(snum-1)))

flipspin(p::Pairing) = typeof(p)(p.level, flip(p.spin))

nlevels(a::TensorBasis) = nlevels(typeof(a))
nlevels(::Type{Pairing{L}}) where L = L

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
