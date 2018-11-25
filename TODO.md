x Fix `@defSub` definition. Since Sub definition can be parameterized, `_makeixmap` has to
  somehow work on-the-fly.
- [IDEA] Redefine (or rather, alias) the type hierarchy as follows:
  ```
  abstract type BasisElem end
  abstract type TensorBasisElem{N} <: BasisElem end
  const TensorBasis{N} = Type{<:TensorBasisElem{N}}
  const Basis{ElemT<:BasisElem} = Type{ElemT}
  ```
  For example, `Pairing` would become
  ```
  basis(b::BasisElem) = typeof(b)
  struct LevelSpin{L} <: TensorBasisElem{2}
      level::Int
      spin::Spin
  end

  nlevels(::Basis{LevelSpin{L}}) where L = L
  nlevels(b::LevelSpin) = nlevel(basis(b))
  index(b::LevelSpin) = ...
  indexbasis(b::Basis{<:LevelSpin}, i, j) = ...
  ```
