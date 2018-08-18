# General
- Move `basistype()`, `rep()`, `reptype` under `ManyBody` root

# Operators
- Factor into multiple files
- Allow `AbstractOperator{N}`s to act on `Product{N}`s

# States
x Interface with Bases
x [TENTATIVE] Algebra
- Add conversion between `ArrayState{}`s of differing `Product{}`s
- [TENTATIVE] Rename `ArrayState{}` to `Tensor{}`

# Bases
- Is there a way to let
  ```
  abstract type Rep{B<:AbstractBasis} end
  ```
that works with Product{1}?
x Add promote_rule between `Product{}`s of `Rep{}`s
- [TENTATIVE] Make distinction between `SimpleBasis{}` and `ComplexBasis{}`
