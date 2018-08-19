# General
- Move `basistype()`, `rep()`, `reptype` under `ManyBody` root
- Refactor code to group together types and individual methods separately.

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
- At least conceptually, this is what we want:
  ```
  abstract type AbstractBasis end
  abstract type ConcreteBasis <: AbstractBasis end
  abstract type AbstractIndex{B<:ConcreteBasis} <: AbstractBasis end

  struct Pairing{L} <: ConcreteBasis end
  struct Slater{B<:ConcreteBasis} <: ConcreteBasis end
  struct Sub{B<:ConcreteBasis, T} <: ConcreteBasis end
  struct Product{B<:NTuple{<:Any, ConcreteBasis}} <: ConcreteBasis end
  struct CartesianIndex{B<:ConcreteBasis, N} <: AbstractIndex{B} end
  struct LinearIndex{B<:ConcreteBasis} <: AbstractIndex{B} end

  abstract type IndexType end
  struct Cartesian{N} <: IndexType end
  struct Linear end

  ### States
  abstract type Mixed{B<:AbstractBasis} <: AbstractState end

  struct Scaled{B<:ConcreteBasis, T} <: Mixed{B} end
  ## Enforce invariants in inner constructors
  struct ArrayState{B<:Bases.AbstractIndex, A<:AbstractArray} <: Mixed{B} end
  const VectorState
  ```
- Create two index types: `LinearIndex` and `CartesianIndex`, similar to what is in `Base`.
