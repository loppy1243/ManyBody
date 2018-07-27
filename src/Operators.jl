@reexport module OperatorsMod
export Operator, SPOperator, MBOperator, tabulate, refop, RaiseOp, LowerOp, RaiseLowerOps, A,
       refop, contract, 

using Combinatorics: permutations, levicivita
using ..States

struct Operator{B, Basis<:SPState, Op}
    op::Op
end
const SPOPerator = Operator
const MBOperator{Basis<:SPState, Op} = Operator{1, Basis, Op}

(op::Operator{B, Basis, <:Function})(sps::Vararg{Basis, 2B}) where {B, Basis} = op.op(sps...)
(op::Operator{B, Basis, <:AbstractArray{T, 2B}})(sps::Vararg{Basis, 2B}) where {B, Basis, T} =
    op.op[map(pnum, sps)...]

tabulate(op::Operator) = tabulate(Float64, op)
tabulate(op::Operator{B, Basis, <:AbstractArray}) where {B, Basis} = op
function tabulate(T::Type, op::Operator{B, Basis, Op}) where {B, Basis, Op}
    states = iter(Basis) |> collect
    statenums = map(pnum, states)
    numstates = length(states)

    mat = Array{T, 2B}(fill(numstates, 2B)...)
    for ixs in Iterators.product(fill(1:numstates, 2B)...)
        ixs = collect(ixs)
        mat[statenums[ixs]...] = op(states[ixs]...)
    end

    Operator{B, Basis, typeof(mat), T}(mat)
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
contract(::Type{R}, a1::RaiseOp{SP}, a2::LowerOp{SP}) =
    a1.state == a2.state ? isocc(R, a1.state) : 0
contract(::Type{R}, a1::LowerOp{SP}, a2::RaiseOp{SP}) =
    a1.state == a2.state ? ~isocc(R, a1.state) : 0
contract(::Type{R}, a::RLOp{SP}) where SP = 0
function contract(::Type{R}, a::RaiseLowerOps) where R<:RefState
    numops = n_ops(a)
    isodd(numops) && return 0

    sum(permuations(1:numops)) do perm
        prod(1:div(numops, 2)) do k
            levicivita(perm)*contract(R, a.ops[perm[2k-1]], a.ops[perm[2k]])
        end
    end
end

(::RLOp{SP})(::MBState{R, SP}, ::MBState{R, SP}) where {R, SP} = 0
function (a::RaiseLowerOps{SP})(X::MBState{R, SP}, Y::MBState{R, SP}) where {R, SP}
    bra_rl_op = refop(X)'
    ket_rl_op - refop(Y)
    
    contract(R, bra_rl_op*a*ket_rl_op)
end

end # module OperatorsMod
