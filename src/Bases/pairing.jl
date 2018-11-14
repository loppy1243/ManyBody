export Pairing, nlevels, flipspin
using ..SpinMod

struct Pairing{Levels} <: TensorBasis{2}
    level::Int
    spin::Spin

    function Pairing{L}(l::Int, s::Spin) where L
        @assert l <= L
        new(l, s)
    end
end

### TensorBasis methods
##############################################################################################
fulldims(B::Type{<:Pairing}) = (nlevels(B), 2)
index(sp::Pairing) = CartesianIndex(p.level, 1+Bool(p.spin))
indexbasis(B::Type{<:Pairing}, level::Int, snum::Int) =
    B(level, Spin(Bool(snum-1)))

### Fermi{<:Pairing} methods
##############################################################################################
nocc(ref::RefStates.Fermi{<:Pairing}) = 2(nlevels(spbasis(f)) - ref.fermilevel)
nunocc(ref::RefStates.Fermi{<:Pairing}) = 2ref.fermilevel
isocc(ref::RefStates.Fermi{<:Pairing}, b) where SPB<:Pairing =
    convert(spbasis(ref)).level <= ref.fermilevel

### Pairing methods
##############################################################################################
SpinMod.spinup(sp::Pairing) = spinup(sp.spin)
flipspin(p::Pairing) = typeof(p)(p.level, flip(p.spin))

nlevels(a::Pairing) = nlevels(typeof(a))
nlevels(::Type{Pairing{L}}) where L = L

Base.show(io::IO, x::Pairing) = print(io, x.level, spinup(x.spin) ? "↑" : "↓")
function Base.show(io::IO, ::MIME"text/plain", x::Pairing)
    L = nlevels(x)
    for i = L:-1:1
        print(io, lpad(i, ndigits(L)), "|  ")
        if i == x.level
            if spindown(x.spin)
                print(io, "↓  ")
            else
                print(io, "o  ")
            end

            if spinup(x.spin)
                println(io, "↑")
            else
                println(io, "o")
            end
        else
            println(io, "o  o")
        end
    end
end
