@reexport module Operators
export Operator, tabulate, refop, RaiseOp, LowerOp, RaiseLowerOps, A, refop, contract,
       set_default_basis!, @default_basis!, @Operator#, @OperatorNBody
using Combinatorics: permutations, levicivita
using Loppy.Util: cartesian_pow
using ..Bases
using ..Bases: Index

#abstract type Operator{B<:Basis, Rep} end
#
#OperatorNBody_sym(N::Int) = Symbol("Operator$(N)Body")
#macro OperatorNBody(N::Int)
#    esc(OperatorNBody_sym(N))
#end
#
#macro def_OperatorNBody(N::Int)
#    _def_OperatorNBody(N, OperatorNBody_sym(N))
#end
#macro def_OperatorNBody(N::Int, sym::Symbol)
#    _def_OperatorNBody(N, sym)
#end
#function _def_OperatorNBody(N::Int, sym::Symbol)
#    ty_sym = esc(sym)
#
#    Basis_sym = esc(gensym("B"))
#
#    sp_syms = [Symbol("sp$i") for i = 1:2N]
#    sp_sym_nums = map(x -> :(index($x)), sp_syms)
#    sp_args = map(x -> :($x::$Basis_sym), sp_syms)
#    sp_sub_args = map(x -> :($x::SubBasis{$Basis_sym}), sp_syms)
#
#    quote
#        struct $ty_sym{B<:Basis, Rep} <: Operator{B, Rep}
#            op::Rep
#        end
#
#        Operators.nbodies(::Type{<:$ty_sym}) = $N
#        (op::$ty_sym{$Basis_sym})($(sp_sub_args...)) where $Basis_sym =
#            op($(map(x -> :($x.state), sp_syms)))
#        (op::$ty_sym{$Basis_sym, <:Function})($(sp_args...)) where $Basis_sym =
#            op.op($(sp_syms...))
#        (op::$ty_sym{$Basis_sym, <:AbstractArray{T, $(2N)}})($(sp_args...)) where
#                {$Basis_sym, T} =
#            op.op[$(sp_sym_nums...)]
#    end
#end
#nbodies() = MethodError(nbodies, ()) |> throw


const SBasis{B<:Basis} = Union{B, <:SubBasis{B}}

struct Operator{N, B<:Basis, Op}
    op::Op
end
@generated (op::Operator{N, B, Op})(args::Vararg{<:SBasis{B}, N2}) where {N, N2, B, Op} =
    :(applyop(Val{$(2N)}, op, $((:(args[$i]) for i=1:N2)...)))
@generated applyop(::Type{Val{N2}}, op::Operator{N, B, <:Function},
                   sps::Vararg{<:SBasis{B}, N2}) where {N, N2, B} =
    :(op.op($((sps[i] !== B && sps[i] <: SubBasis ? :(sps[$i].state) : :(sps[$i])
               for i = 1:N2)...)))
@generated applyop(::Type{Val{N2}}, op::Operator{N, B, <:AbstractArray{T, N2}},
                   sps::Vararg{Union{B, SB}, N2}) where {N, N2, B, T, SB<:SubBasis{B}} =
    :(op.op[$((sps[i] !== B && <: SubBasis ? :(index(sps[$i].state)) : :(index(sps[$i]))
               for i = 1:N2)...)])


#for i = 1:3
#    @eval @def_OperatorNBody $i
#    @eval export $(OperatorNBody_sym(i))
#end

nbodies(::Type{<:Operator{N}}) where N = N
nbodies(op::Operator) = nbodies(typeof(op))

tabulate(op) = tabulate(Complex64, op)
tabulate(op::Operator{N, <:Index, <:AbstractArray}) where N = op
function tabulate(::Type{T}, op::Operator{N, B, Op}) where {T, N, B, Op}
    mat = Array{T}(fill(dim(B), 2N)...)
    for ss in cartesian_pow(B, 2N)
        mat[CartesianIndex(map(index, ss))] = op(ss...)
    end

    Operator{N, Index{B}, typeof(mat)}(mat)
end

struct Raised{T}; val::T end
↑(x) = Raised(x)

struct RaiseOp{SP<:SPBasis}; state::SP end
struct LowerOp{SP<:SPBasis}; state::SP end
RaiseOp(s::SubBasis{SP}) where SP = RaiseOp{SP}(s.state)
LowerOp(s::SubBasis{SP}) where SP = LowerOp{SP}(s.state)
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

const SPSubBasis{SP<:SPBasis} = Union{SP, <:SubBasis{SP}}
const MBSubBasis{MB<:MBBasis} = Union{MB, <:SubBasis{MB}}

A(p::SPSubBasis) = LowerOp(p)
A(p::Bra{<:SPSubBasis}) = RaiseOp(p.state)

#A(::Type{SP}, sps::NTuple{N, Union{T, Bra{T}, Raised{T}}}) where {SP, N, T<:SPSubBasis{SP}} =
#    A(sps)
@generated function A(sps::NTuple{N, Union{T, Bra{T}, Raised{T}}}) where {SP, T<:SPSubBasis{SP}, N}
    args = map(enumerate(sps.parameters)) do x
        i, T = x
        if T <: Bra
            :(RaiseOp(sps[$i].state))
        elseif T <: Raised
            :(RaiseOp(sps[$i].val))
        else
            :(LowerOp(sps[$i]))
        end
    end

    :(RaiseLowerOps{SP}([$(args...)]))
end
#@generated function A(::Type{SP},
#                      sps::NTuple{N, Union{Int, Index{SP}, Raised{Int},
#                                           Raised{Index{SP}}}}) where {SP, N}
#    args = map(enumerate(sps.parameters)) do x
#        i, T = x
#
#        if T <: Index
#            :(LowerOp{Index{SP}}(sps[$i]))
#        elseif T == Int
#            :(LowerOp{Index{SP}}(sps[$i]))
#        elseif T == Raised{Int}
#            :(RaiseOp{Index{SP}}(sps[$i].val))
#        else
#            :(RaiseOp{Index{SP}}(sps[$i]))
#        end
#    end
#
#    :(RaiseLowerOps{Index{SP}}([$(args...)]))
#end

n_ops(::RLOp) = 1
n_ops(a::RaiseLowerOps) = length(a.ops)

refop(s::SubBasis{B}) where B<:MBBasis = refop(s.state)
function refop(s::Bases.PartHole{R}) where {SP, R<:RefState{SP}}
    a = RaiseLowerOps{SP}(map(x -> RaiseOp{SP}(indexp(R, x)), find(s.parts)))
    b = RaiseLowerOps{SP}(map(x -> LowerOp{SP}(indexh(R, x)), find(s.holes)))
    a*b
end

contract(::Type, a1::OP, a2::OP) where OP<:RLOp = 0
contract(::Type{Bases.PartHole{R}}, a1::RaiseOp{SP}, a2::LowerOp{SP}) where
        {SP, R<:RefState{SP}} =
    a1.state == a2.state ? isocc(R, a1.state) : 0
contract(::Type{Bases.PartHole{R}}, a1::LowerOp{SP}, a2::RaiseOp{SP}) where
        {SP, R<:RefState{SP}} =
    a1.state == a2.state ? ~isocc(R, a1.state) : 0
contract(::Type, a::RLOp) = 0

function contract(MB::Type, a::RaiseLowerOps)
    numops = n_ops(a)
    isodd(numops) && return 0

    sum(permutations(1:numops)) do perm
        prod(1:div(numops, 2)) do k
            levicivita(perm)*contract(MB, ops[perm[2k-1]], ops[perm[2k]])
        end
    end
end

for op_ty in [RaiseOp, LowerOp]
    @eval (::$op_ty)(::MB, ::MB) where MB<:MBBasis = 0
end

function (a::RaiseLowerOps)(X::MBSubBasis{MB}, Y::MBSubBasis{MB}) where MB<:MBBasis
    bra_rl_op = refop(X)'
    ket_rl_op = refop(Y)
    
    contract(MB, bra_rl_op*a*ket_rl_op)
end

DEFAULT_BASIS = nothing
DEFAULT_REFSTATE = nothing

set_default_basis!(::Type{B}) where B<:Basis = global DEFAULT_BASIS = B

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
        Operator{$nbodies, $(esc(basis)), typeof(x)}(x)
    end
end

end # module Operators
