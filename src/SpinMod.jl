@reexport module SpinMod
export Spin, SPINUP, SPINDOWN, ↑, ↓, SPINS, flip, spinup, spindown

@enum Spin::Bool SPINUP=true SPINDOWN=false
const ↑, ↓ = SPINUP, SPINDOWN

const SPINS = [SPINDOWN, SPINUP]

flip(s::Spin) = Spin(~Bool(s))
spinup(s::Spin) = Bool(s)
spindown(s) = ~spinup(s)

#Base.show(io::IO, x::Spin) = println(io, spinup(x) ? '↑' : '↓')

end # module SpinMod
