abstract type AbstractState end

module Bases

abstract type AbstractBasis <: AbstractState
abstract type ConcreteBasis <: AbstractBasis
abstract type AbstractIndex{B<:ConcreteBasis} <: AbstractBasis

#module GenMethods
#    export Generation
#
#    abstract type Generation end
#    struct Computed <: Generation end
#    struct Cached <: Generation end
#    struct Provided <: Generation end
#end
#using .GenMethod

struct Pairing{L} <: ConcreteBasis
    level::Int
    spin::Spin

    function Pairing{L}(level, spin::Spin) where L
        level::Int = level
        @assert level <= L
        new(level, spin)
    end
end
#Generation(::Type{<:Pairing}) = GenMethods.Cached()

level(b::Pairing) = b.level
spin(b::Pairing) = b.spin
flipspin(b::Pairing) = typeof(b)(level(b), flipspin(spin(b)))
nlevels(::Type{Pairing{L}}) where L = L

dim(B::Type{<:Pairing}) = 2nlevels(B)
function indexbasis(B::Type{<:Pairing}, i::Int)
end

end # module Bases



module States

end # module States
