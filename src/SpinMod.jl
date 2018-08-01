@reexport module SpinMod
export Spin, SPINUP, SPINDOWN, SPINS, flip, spinup, spindown

@enum Spin::Bool SPINUP=true SPINDOWN=false

const SPINS = [SPINDOWN, SPINUP]

flip(s::Spin) = Spin(~Bool(s))
spinup(s::Spin) = Bool(s)
spindown(s::Spin) = ~Bool(s)

end # module SpinMod
