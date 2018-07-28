@reexport module Operators
export Operator, SPOperator, MBOperator, tabulate, refop, RaiseOp, LowerOp, RaiseLowerOps, A,
       refop, contract

using Combinatorics: permutations, levicivita
using ..States

abstract type Operator{Basis<:SPState, Rep} end

macro def_OperatorNBody(N::Int)
    _def_OperatorNBody(N, Symbol("Operator$(N)Body"))
end
macro def_OperatorNBody(N::Int, sym::Symbol)
    _def_OperatorNBody(N, sym)
end
function _def_OperatorNBody(N::Int, sym::Symbol)
    ty_sym = esc(sym)

    Basis_sym = esc(gensym("Basis"))

    sp_syms = [Symbol("sp$i") for i = 1:2N]
    sp_sym_nums = map(x -> :(x.num), sp_syms)
    sp_args = map(x -> :($x::$Basis_sym), sp_syms)

    quote
        struct $ty_sym{Basis<:SPState, Rep} <: Operator{Basis, Rep}
            op::Rep
        end

        Operators.nbodies(::Type{<:$ty_sym}) = $N
        (op::$ty_sym{$Basis_sym, <:Function})($(sp_args...)) where $Basis_sym =
            op.op($(sp_syms...))
        (op::$ty_sym{$Basis_sym, <:Function})($(sp_args...)) where $Basis_sym<:IntState =
            op.op($(sp_sym_nums...))
        (op::$ty_sym{$Basis_sym, <:AbstractArray{T, $(2N)}})($(sp_args...)) where
                {$Basis_sym<:IntState, T} =
            op.op[$(sp_sym_nums...)]
    end
end
nbodies() = MethodError(nbodies, ()) |> throw

for i = 1:3; @eval @def_OperatorNBody $i; end
@def_OperatorNBody 1 MBOperator

nbodies(op::Operator) = nbodies(typeof(op))

tabulate(op::Operator) = tabulate(Float64, op)
tabulate(op::Operator{Basis, <:AbstractArray}) where Basis<:IntState = op
function tabulate(T::Type, op::Operator{Basis, Op}) where {Basis, Op}
    states = iter(Basis) |> collect
    statenums = map(snum, states)
    numstates = length(states)

    mat = Array{T, 2B}(fill(numstates, 2B)...)
    for ixs in Iterators.product(fill(1:numstates, 2B)...)
        ixs = collect(ixs)
        mat[statenums[ixs]...] = op(states[ixs]...)
    end

    Operator{IntState{Basis}, typeof(mat)}(mat)
end

struct RaiseOp{SP<:SPState}; state::SP end
struct LowerOp{SP<:SPState}; state::SP end
struct RaiseLowerOps{SP<:SPState}
    ops::Vector{Union{RaiseOp{SP}, LowerOp{SP}}}
end
RaiseLowerOps(sps::Vararg{Union{SP, Bra{SP}}}) where SP =
    RaiseLowerOps{SP}([sp isa Bra ? RaiseOp(snum(sp.state)) : LowerOp(snum(sp)) for sp in sps])
RaiseLowerOps(itr) =
    RaiseLowerOps{typeof(first(itr))}([sp isa Bra ? RaiseOp(snum(sp.state)) : LowerOp(snum(sp))
                                      for sp in itr])

const RLOp{SP} = Union{RaiseOp{SP}, LowerOp{SP}}
Base.:*(op1::RLOp{SP}, op2::RLOp{SP}) where SP = RaiseLowerOps{SP}([op1, op2])
Base.:*(op1::RLOp{SP}, op2::RaiseLowerOps{SP}) where SP = RaiseLowerOps{SP}([op1; op2.ops])
Base.:*(op1::RaiseLowerOps{SP}, op2::RLOp{SP}) where SP = RaiseLowerOps{SP}([op2.ops; op2])
Base.:*(op1::RaiseLowerOps{SP}, op2::RaiseLowerOps{SP}) where SP =
    RaiseLowerOps{SP}([op1.ops; op2.ops])
Base.ctranspose(op::RaiseOp{SP}) where SP = LowerOp{SP}(op.state)
Base.ctranspose(op::LowerOp{SP}) where SP = RaiseOp{SP}(op.state)
Base.ctranspose(op::RaiseLowerOps{SP}) where SP =
    RaiseLowerOps{SP}(reverse(map(ctranspose, op.ops)))
## Don't know if we want this.
#Base.getindex(op::RaiseLowerOps, i) = typeof(op)([op.ops[i]])

A(p::SPState) = LowerOp{typeof(p)}(p)
A(p::Bra) = RaiseOp{typeof(p.state)}(p.state)
A(x) = RaiseLowerOps(x)

n_ops(::RLOp) = 1
n_ops(a::RaiseLowerOps) = length(a.ops)

refop(s::MBState{R, SP}) where {R, SP} =
    A(map(x -> Bra(IntState{SP}(x)), find(s.parts))) * A(map(IntState{SP}, find(s.holes)))

contract(::Type{<:RefState}, a1::OP, a2::OP) where OP<:RLOp = 0
contract(::Type{Vacuum{SP}}, a1::RaiseOp{SP}, a2::LowerOp{SP}) where SP = 0
contract(::Type{Vacuum{SP}}, a1::LowerOp{SP}, a2::RaiseOp{SP}) where SP =
    a1.state == a2.state
contract(::Type{R}, a1::RaiseOp{SP}, a2::LowerOp{SP}) where {R, SP}=
    a1.state == a2.state ? isocc(R, a1.state) : 0
contract(::Type{R}, a1::LowerOp{SP}, a2::RaiseOp{SP}) where {R, SP} =
    a1.state == a2.state ? ~isocc(R, a1.state) : 0
contract(::Type{R}, a::RLOp{SP}) where {R, SP} = 0
function contract(::Type{R}, a::RaiseLowerOps) where R<:RefState
    numops = n_ops(a)
    isodd(numops) && return 0

    sum(permuations(1:numops)) do perm
        prod(1:div(numops, 2)) do k
            levicivita(perm)*contract(R, a.ops[perm[2k-1]], a.ops[perm[2k]])
        end
    end
end

for op_ty in [RaiseOp, LowerOp]
    @eval (::$op_ty{SP})(::MBState{R, SP}, ::MBState{R, SP}) where {R, SP} = 0
end

function (a::RaiseLowerOps{SP})(X::MBState{R, SP}, Y::MBState{R, SP}) where {R, SP}
    bra_rl_op = refop(X)'
    ket_rl_op - refop(Y)
    
    contract(R, bra_rl_op*a*ket_rl_op)
end

end # module Operators
