@reexport module Operators
export Operator, SPOperator, MBOperator, tabulate, refop, RaiseOp, LowerOp, RaiseLowerOps, A,
       refop, contract, set_default_basis!, set_default_refstate!, @default_basis!,
       @default_refstate!, @Operator, @MBOperator, @OperatorNBody, OperatorNBody_sym

using Combinatorics: permutations, levicivita
using Loppy.Util: cartesian_pow
using ..Bases

abstract type Operator{B<:Basis, Rep} end

OperatorNBody_sym(N::Int) = Symbol("Operator$(N)Body")
macro OperatorNBody(N::Int)
    esc(OperatorNBody_sym(N))
end

macro def_OperatorNBody(N::Int)
    _def_OperatorNBody(N, OperatorNBody_sym(N))
end
macro def_OperatorNBody(N::Int, sym::Symbol)
    _def_OperatorNBody(N, sym)
end
function _def_OperatorNBody(N::Int, sym::Symbol)
    ty_sym = esc(sym)

    Basis_sym = esc(gensym("B"))

    sp_syms = [Symbol("sp$i") for i = 1:2N]
    sp_sym_nums = map(x -> :(index($x)), sp_syms)
    sp_args = map(x -> :($x::$Basis_sym), sp_syms)

    quote
        struct $ty_sym{B<:Basis, Rep} <: Operator{B, Rep}
            op::Rep
        end

        Operators.nbodies(::Type{<:$ty_sym}) = $N
        (op::$ty_sym{$Basis_sym, <:Function})($(sp_args...)) where $Basis_sym =
            op.op($(sp_syms...))
        (op::$ty_sym{$Basis_sym, <:Function})($(sp_args...)) where $Basis_sym<:Index =
            op.op($(sp_sym_nums...))
        (op::$ty_sym{$Basis_sym, <:AbstractArray{T, $(2N)}})($(sp_args...)) where
                {$Basis_sym, T} =
            op.op[$(sp_sym_nums...)]
    end
end
nbodies() = MethodError(nbodies, ()) |> throw

for i = 1:3
    @eval @def_OperatorNBody $i
    @eval export $(OperatorNBody_sym(i))
end
@def_OperatorNBody 1 MBOperator

nbodies(op::Operator) = nbodies(typeof(op))

tabulate(op::Operator) = tabulate(Complex64, op)
tabulate(op::Operator{B, <:AbstractArray}) where B<:Index = op
function tabulate(T::Type, op::Operator{B, Op}) where {B, Op}
    n = 2nbodies(typeof(op))
    mat = Array{T}(fill(dim(B), n)...)
    for ixs in cartesian_pow(indices(B), n)
        mat[CartesianIndex(ixs)] = op(map(i -> B[i], ixs)...)
    end

    Operator{Index{B}, typeof(mat)}(mat)
end

struct Raised{T}; val::T end
↑(x) = Raised(x)

struct RaiseOp{SP<:SPBasis}; state::SP end
struct LowerOp{SP<:SPBasis}; state::SP end
const RLOp{SP} = Union{RaiseOp{SP}, LowerOp{SP}}

struct RaiseLowerOps{SP<:SPBasis}
    ops::Vector{RLOp{SP}}
end

Base.convert(::Type{RaiseLowerOps{SP}}, a::RLOp{SP}) where SP = RaiseLowerOps{SP}([a])
Base.convert(::Type{RaiseLowerOps}, a::RLOp{SP}) where SP = RaiseLowerOps{SP}([a])

Base.:*(op1::RLOp{SP}, op2::RLOp{SP}) where SP = RaiseLowerOps{SP}([op1, op2])
Base.:*(op1::RLOp{SP}, op2::RaiseLowerOps{SP}) where SP = RaiseLowerOps{SP}([op1; op2.ops])
Base.:*(op1::RaiseLowerOps{SP}, op2::RLOp{SP}) where SP = RaiseLowerOps{SP}([op2.ops; op2])
Base.:*(op1::RaiseLowerOps{SP}, op2::RaiseLowerOps{SP}) where SP =
    RaiseLowerOps{SP}([op1.ops; op2.ops])
Base.ctranspose(op::RaiseOp{SP}) where SP = LowerOp{SP}(op.state)
Base.ctranspose(op::LowerOp{SP}) where SP = RaiseOp{SP}(op.state)
Base.ctranspose(op::RaiseLowerOps{SP}) where SP =
    RaiseLowerOps{SP}(reverse(map(ctranspose, op.ops)))

A(sps...) = A(sps)
A(T::Type, sps...) = A(T, sps)

A(p::SPBasis) = LowerOp{typeof(p)}(p)
A(p::Bra{<:SPBasis}) = RaiseOp{typeof(p.state)}(p.state)
A(::Type{SP}, p::SP) where SP<:SPBasis = A(p)
A(::Type{SP}, p::Bra{SP}) where SP<:SPBasis = A(p)

A(::Type{SP}, sps::Tuple{Vararg{Union{SP, Bra{SP}, Raised{SP}}, N}}) where {SP, N} = A(sps)
@generated function A(sps::Tuple{Vararg{Union{SP, Bra{SP}, Raised{SP}}, N}}) where {SP, N}
    args = map(enumerate(sps.parameters)) do x
        i, T = x
        if T <: Bra
            :(RaiseOp{SP}(sps[$i].state))
        elseif T <: Raised
            :(RaiseOp{SP}(sps[$i].val))
        else
            :(LowerOp{SP}(sps[$i]))
        end
    end

    :(RaiseLowerOps{SP}([$(args...)]))
end
@generated function A(::Type{SP},
                      sps::Tuple{Vararg{Union{Int, Index{SP}, Raised{Int},
                                              Raised{Index{SP}}}, N}}) where {SP, N}
    args = map(enumerate(sps.parameters)) do x
        i, T = x

        if T <: Index
            :(LowerOp{Index{SP}}(sps[$i]))
        elseif T == Int
            :(LowerOp{Index{SP}}(sps[$i]))
        elseif T == Raised{Int}
            :(RaiseOp{Index{SP}}(sps[$i].val))
        else
            :(RaiseOp{Index{SP}}(sps[$i]))
        end
    end

    :(RaiseLowerOps{Index{SP}}([$(args...)]))
end

n_ops(::RLOp) = 1
n_ops(a::RaiseLowerOps) = length(a.ops)

function refop(s::MBBasis{R, SP}) where {R, SP}
    a = RaiseLowerOps{SP}(map(x -> RaiseOp{SP}(indexp(R, x)), find(s.parts)))
    b = RaiseLowerOps{SP}(map(x -> LowerOp{SP}(indexh(R, x)), find(s.holes)))
    a*b
end

contract(::Type, a1::OP, a2::OP) where OP<:RLOp = 0
contract(R::Type, a1::RaiseOp{SP}, a2::LowerOp{SP}) where SP =
    a1.state == a2.state ? isocc(R, a1.state) : 0
contract(R::Type, a1::LowerOp{SP}, a2::RaiseOp{SP}) where SP =
    a1.state == a2.state ? ~isocc(R, a1.state) : 0
contract(::Type, a::RLOp) = 0
function contract(R::Type, a::RaiseLowerOps)
    numops = n_ops(a)
    isodd(numops) && return 0

    sum(permutations(1:numops)) do perm
        prod(1:div(numops, 2)) do k
            levicivita(perm)*contract(R, a.ops[perm[2k-1]], a.ops[perm[2k]])
        end
    end
end

for op_ty in [RaiseOp, LowerOp]
    @eval (::$op_ty{SP})(::MBBasis{R, SP}, ::MBBasis{R, SP}) where {R, SP} = 0
end

function (a::RaiseLowerOps{SP})(X::MBBasis{R, SP}, Y::MBBasis{R, SP}) where {R, SP}
    bra_rl_op = refop(X)'
    ket_rl_op = refop(Y)
    
    contract(R, bra_rl_op*a*ket_rl_op)
end

DEFAULT_BASIS = nothing
DEFAULT_REFSTATE = nothing

set_default_basis!(::Type{B}) where B<:Basis = global DEFAULT_BASIS = B
set_default_refstate!(::Type{R}) where R<:RefState = global DEFAULT_REFSTATE = R

macro default_basis!(sym::Symbol)
    :(set_default_basis!($(esc(sym))))
end
macro default_basis!(expr::Expr)
    if expr.head == :curly || expr.head == :.
        return :(set_default_basis!($(esc(expr))))
    end
    @assert expr.head == :(=) && (expr.args[1] isa Symbol || expr.args[1].head != :call) #=
         =# || expr.head == :const

    is_const = expr.head == :const
    name = is_const ? expr.args[1].args[1] : expr.args[1]

    expr = is_const ? Expr(:const, expr) : expr

    quote
        $(esc(expr))
        set_default_basis!($(esc(name)))
    end
end

macro default_refstate!(sym::Symbol)
    :(set_default_refstate!($(esc(sym))))
end
macro default_refstate!(expr::Expr)
    if expr.head == :curly || expr.head == :.
        return :(set_default_refstate!($(esc(expr))))
    end
    @assert expr.head == :(=) && (expr.args[1] isa Symbol || expr.args[1].head != :call) #=
         =# || expr.head == :const

    is_const = expr.head == :const
    name = is_const ? expr.args[1].args[1] : expr.args[1]

    expr = is_const ? Expr(:const, expr) : expr
    quote
        $(esc(expr))
        set_default_refstate!($(esc(name)))
    end
end

macro Operator(expr)
    :(@Operator($(esc(expr)), DEFAULT_BASIS))
end
macro Operator(expr::Expr, basis)
    @assert expr.head == :(->) && expr.args[1].head == :tuple
    nbodies2 = length(expr.args[1].args)
    @assert iseven(nbodies2)
    nbodies = div(nbodies2, 2)

    quote
        x = $(esc(expr))
        @OperatorNBody($nbodies){$(esc(basis)), typeof(x)}(x)
    end
end

macro MBOperator(expr)
    :(@MBOperator($(esc(expr)), DEFAULT_REFSTATE, DEFAULT_BASIS))
end
macro MBOperator(expr::Expr, refstate, basis)
    @assert expr.head == :(->) && expr.args[1].head == :tuple
    @assert length(expr.args[1].args) == 2

    quote
        x = $(esc(expr))
        MBOperator{MBBasis{$(esc(refstate)), $(esc(basis))}, typeof(x)}(x)
    end
end

end # module Operators
