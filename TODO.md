# General
- Move `basistype()`, `rep()`, `reptype` under `ManyBody` root
  - [NOTE] Moved `basistype()`
- Refactor code to group together types and individual methods separately.
- Move `RefStates` out of `Bases`
- Update tests to match new scheme
- [IDEA] Move `dim()`, `rank()` up to `ManyBody`?
- [IDEA] Define within a given module only that which belongs to the module, and define
  interactions or outside methods within the greater module.

# Operators
- Factor into multiple files
- Allow `AbstractOperator{N}`s to act on `Product{N}`s

# States
x Interface with Bases
x [IDEA] Algebra
x [IGNORE] Add conversion between `ArrayState{}`s of differing `Product{}`s
  - [NOTE] Do *NOT* do this, instead only worry about exact types. Conversions can happen
    upon action.
- [IDEA] Rename `ArrayState{}` to `Tensor{}`
- Add methods for `similar()` for `ArrayState{}`

# Bases
- Is there a way to let
  ```
  abstract type Rep{B<:AbstractBasis} end
  ```
that works with Product{1}?
x Add promote_rule between `Product{}`s of `Rep{}`s
x Make distinction between `SimpleBasis{}` and `ComplexBasis{}`
  - [NOTE] See below item for actual scheme implemented
x At least conceptually, this is what we want:
  ```
  abstract type AbstractBasis end
  abstract type ConcreteBasis <: AbstractBasis end
  abstract type LeafBasis <: AbstractBasis end
  abstract type AbstractIndex{B<:ConcreteBasis} <: LeafBasis end

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
x Create two index types: `LinearIndex` and `CartesianIndex`, similar to what is in `Base`.
- Add way to select `IndexType(::Type{<:Sub})` in `@defSub`
  - Define `innerdims()`
