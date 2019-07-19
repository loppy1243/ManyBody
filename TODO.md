* [ ] [IDEA] Redefine (or rather, alias) the type hierarchy as follows:
```julia
abstract type BasisElem end
abstract type TensorBasisElem{N} <: BasisElem end
abstract type AbstractBasis end
struct Basis{ElemT<:BasisElem} <: AbstractBasis end
const TensorBasis{N, ElemT<:TensorBasisElem{N}} = Basis{ElemT}
```
For example, `Pairing` would become
```julia
basis(b::BasisElem) = typeof(b)
struct LevelSpin{L} <: TensorBasisElem{2}
    level::Int
    spin::Spin
end

nlevels(::Basis{LevelSpin{L}}) where L = L
nlevels(b::LevelSpin) = nlevels(basis(b))
index(b::LevelSpin) = ...
indexbasis(b::Basis{<:LevelSpin}, i, j) = ...
```
* [x] Remove `tabulate` -- it doesn't do anything broadcasting can't.
* [x] Stop punning on the basis type as an iterator.
* [-] Remove `SubBasis`, replace with i.e. `FilteredBasis`, `MappedBasis`, etc.
* [ ] Remove `convert(::Type{<:Array}), ::AbstractBasis)` methods, replace with constructors
  methods.
* [ ] Make `DualBasis` sufficiently general; in particular, consider the rank > 2 case.
* [ ] [IDEA] Concerning `SubBasis`: If we adopt the idea above with seperated `Basis` and
  `BasisElem`, then we can pass around `BasisElem` no problem while defining `Sub` off of
  `Basis`. Something like:
```julia
struct SubBasis{B<:Basis, A<:Array{Int}} <: AbstractBasis
    idxmap::A

    function SubBasis{B}(idxmap::Vector{Int}) where B<:AbstractBasis
        new{B, Vector{Int}}(idxmap)
    end
    function SubBasis{B}(idxmap::Array{Int, M}) where {M, B<:TensorBasis}
        @assert M <= rank(B)
        new{B, Array{Int, M}}(idxmap)
    end
end

## Compare to similar definition in src/Bases/main.jl
paired(A, L) = SubBasis{Slater{Pairing{L}, 2}}() do s
    cnt = 0
    for p in eachocc(s)
        cnt += flipspin(p) in s || return false
    end
    cnt == A
end
```

# Done
* [x] Fix `@defSub` definition. Since Sub definition can be parameterized, `_makeixmap` has to
  somehow work on-the-fly.
