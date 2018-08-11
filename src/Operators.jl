@reexport module Operators
export Operator, tabulate, refop, RaiseOp, LowerOp, RaiseLowerOps, A, refop, contract,
       normord, apply_normord_rl, @Operator

using Base.Cartesian
using Combinatorics: levicivita
using JuliaUtil: cartesian_pow
using ..Bases
using ..States

abstract type AbstractOperator{N, B<:AbstractBasis, T} end
struct FunctionOperator{N, B<:AbstractBasis, T, F<:Function} <: AbstractOperator{N, B, T}
    rep::F
end
struct ArrayOperator{N, B<:AbstractBasis, T, A<:AbstractArray{T, N}} #=
    =# <: AbstractOperator{N, B, T}
    rep::A
end
FunctionOperator{N, B, T}(f::Function) where {N, B<:AbstracatBasis, T} =
    FunctionOperator{N, B, T, typeof(f)}(f)
ArrayOperator{N, B<:AbstractBasis}(a::AbstractArray) where {N, B<:AbstractBasis} =
    ArrayOperator{N, B, eltype(a), typeof(a)}(a)

rep(op::AbstractOperator) = op.rep
Base.eltype(::Type{<:AbstractOperator{<:Any, <:Any, T}}) where T = T

(op::AbstractOperator{N})(::Type{T}, args...)
(op::AbstractOperator{N})(::Type{T}, args...) where {N, T} = op(Array{T, N}, args...)
@generated function (op::AbstractOperator{N, B})(::Type{A}, arg::Vararg{Any, N}) where
                    {N, B<:AbstractBasis, A<:AbstractArray{T, N}}
    X_syms = [Symbol("X_$i") for i = 1:N]
    Y_syms = [Symbol("Y_$i") for i = 1:N]
    quote
        ret = A(undef, $(fill(dim(B), N)...))
        @nloops $N X (_ -> B) begin
            @nloops $N Y (_ -> B) begin
                @nref($N, ret, d -> index(X_d)) +=
                    op[$(X_syms...), $(Y_syms...)] * @ncall($N, +, d -> Y_d'args[d])
            end
        end

        ret
    end
end

@generated function Base.getindex(op::AbstractOperator{N, B}, args::Vararg{B, N2}) where
                                 {N, N2, B<:AbstractBasis}
    @assert N2 == 2N
    :(@ncall($N2, matrixelem, op, i -> args[i]))
end

@generated matrixelem(op::ArrayOperator{N, B}, args::Vararg{B, N2}) where
                     {N, N2, B<:AbstractBasis} =
    :(@nref($N2, op.arr, i -> index(args[i])))
@generated matrixelem(op::FunctionOperator{N, B, T}, args::Vararg{B, N2}) where
                     {N, N2, T, B<:AbstractBasis} =
    :(convert(T, @ncall($N2, op.func, i -> args[i])))
@generated matrixelem(op::GroupedOperator{N, B}, args::Vararg{B, N2}) where
                     {N, N2, B<:AbstractBasis} = quote
    i, j = @ncall($N2, groupindex, op.map, i -> args[i])
    op.mat[i, j]
end
@generated groupindex(f::Function, args::Vararg{<:Any, N}) where N =
    :(@ncall($N, f, i -> args[i]))
@generated groupindex(a::AbstractArray{<:Any, N}, args::Vararg{<:Any, N}) where N =
    :(@nref($N, a, i -> index(args[i])))

nbodies(::Type{<:AbstractOperator{N}}) where N = N
nbodies(op::AbstractOperator) = nbodies(typeof(op))

tabulate(op) = tabulate(ComplexF64, op)
tabulate(::Type, x::Number) = x
tabulate(::Type, x::AbstractArray) = x
tabulate(::Type{T}, op::ArrayOperator{N, B}) where {T, N, B<:Bases.Index} =
    ArrayOperator{N, B}(convert(AbstractArray{T}, op.arr))
tabulate(::Type{T}, op::ArrayOperator{N, B}) where {T, N, B<:AbstractBasis} =
    ArrayOperator{N, Bases.Index{B}}(convert(AbstractArray{T}, op.arr))
tabulate(::Type{T}, op::GroupedOperator{N, B}) where {T, N, B<:Bases.Index} =
    GroupedOperator{N, B}(convert(AbstractMatrix{T}, op.mat), op.map)
tabulate(::Type{T}, op::GroupedOperator{N, B}) where {T, N, B<:AbstractBasis} =
    GroupedOperator{N, Bases.Index{B}}(convert(AbstractMatrix{T}, op.mat), op.map)
@generated function tabulate(::Type{T}, op::AbstractOperator{N, B}) where
           {T, N, B<:AbstractBasis}
    quote
        arr = @ncall($(2N), Array{T}, undef, _ -> dim(B))
        @nloops $(2N) s (_ -> B) begin
            @nref($(2N), arr, i -> index(s_i)) = @nref($(2N), op, s)
        end

        ArrayOperator{N, Bases.Index{B}}(mat)
    end
end

@generated group(f, op::ArrayOperator{N, B}) where {N, B<:AbstractBasis} = quote
    mat = similar(op.arr, dim(B), dim(B))
    @nloops $(2N) b (_ -> B) begin
        mat[CartesianIndex(@ncall($(2N), groupindex, f, b))] = @nref($(2N), op, b)
    end

    GroupedOperator{N, B}(mat, f)
end
@generated ungroup(op::GroupedOperator{N, B}) where {N, B<:AbstractBasis} = quote
    arr = @ncall($(2N), similar, op.mat, _ -> dim(B))
    @nloops $(2N) b (_ -> B) begin
        @nref($(2N), arr, i -> b_i) = mat[CartesianIndex(@ncall($(2N), groupindex, op.map, b))]
    end

    ArrayOperator{N, B}(arr)
end

struct Raised{T}; val::T end
↑(x) = Raised(x)

struct RaiseOp{B<:AbstractBasis}; state::B end
struct LowerOp{B<:AbstractBasis}; state::B end
RaiseOp(s::Bases.Sub{B}) where B = RaiseOp{B}(s.state)
LowerOp(s::Bases.Sub{B}) where B = LowerOp{B}(s.state)
const RLOp{B} = Union{RaiseOp{B}, LowerOp{B}}

struct RaiseLowerOps{B<:AbstractBasis}
    ops::Vector{RLOp{B}}
end

Base.:(==)(::RLOp, RLOp) = false
Base.:(==)(a::R, b::R) where R<:RLOp = a.state == b.state
Base.:(==)(::RaiseLowerOps, ::RaiseLowerOps) = false
Base.:(==)(a::R, b::R) where R<:RaiseLowerOps = a.ops == b.ops

Base.convert(::Type{RaiseLowerOps{B}}, a::RLOp{B}) where B<:AbstractBasis =
    RaiseLowerOps{B}([a])
Base.convert(::Type{RaiseLowerOps}, a::RLOp{B}) where B<:AbstractBasis = RaiseLowerOps{B}([a])

Base.:*(op1::RLOp{B}, op2::RLOp{B}) where B<:AbstractBasis = RaiseLowerOps{B}([op1, op2])
Base.:*(op1::RLOp{B}, op2::RaiseLowerOps{B}) where B<:AbstractBasis =
    RaiseLowerOps{B}([op1; op2.ops])
Base.:*(op1::RaiseLowerOps{B}, op2::RLOp{B}) where B<:AbstractBasis =
    RaiseLowerOps{B}([op2.ops; op2])
Base.:*(op1::RaiseLowerOps{B}, op2::RaiseLowerOps{B}) where B<:AbstractBasis =
    RaiseLowerOps{B}([op1.ops; op2.ops])
Base.adjoint(op::RaiseOp{B}) where B<:AbstractBasis = LowerOp{B}(op.state)
Base.adjoint(op::LowerOp{B}) where B<:AbstractBasis = RaiseOp{B}(op.state)
Base.adjoint(op::RaiseLowerOps{B}) where B<:AbstractBasis =
    RaiseLowerOps{B}(reverse(map(adjoint, op.ops)))

A(sps...) = A(sps)
A(T::Type, sps...) = A(T, sps)

A(p::AbstractBasis) = LowerOp(p)
A(p::Bra{<:AbstractBasis}) = RaiseOp(p.state)

@generated function A(sps::NTuple{N, Union{T, Bra{T}, Raised{T}}}) where
                     {B<:AbstractBasis, T<:Bases.MaybeSub{B}, N}
    args = map(enumerate(sps.parameters)) do x
        i, U = x
        if U <: Bra
            :(RaiseOp(sps[$i].state))
        elseif U <: Raised
            :(RaiseOp(sps[$i].val))
        else
            :(LowerOp(sps[$i]))
        end
    end

    :(RaiseLowerOps{B}([$(args...)]))
end

n_ops(::RLOp) = 1
n_ops(a::RaiseLowerOps) = length(a.ops)

refop(s::Bases.Sub{B}) where B<:AbstractBasis = refop(s.state)
function refop(s::Bases.Slater{B}) where B<:AbstractBasis
    a = RaiseLowerOps{B}(map(x -> RaiseOp{B}(indexp(R, x)), find(s.parts)))
    b = RaiseLowerOps{B}(map(x -> LowerOp{B}(indexh(R, x)), find(s.holes)))
    a*b
end

normord(a::RaiseLowerOps{B}) where B<:AbstractBasis = normord(RefStates.Vacuum{B}, a)
function normord(::Type{R}, a::RaiseLowerOps{B}) where {B<:AbstractBasis, R<:RefState{B}}
    function comp(a, b)
        if a[2] && b[2]
            index(a[1]) >= index(b[1])
        elseif !a[2] && !b[2]
            index(a[1]) <= index(b[1])
        else
            a[2]
        end
    end

    xs = map(a.ops) do b
        (b.state, if b isa RaiseOp
            ~isocc(R, b.state)
        else
            isocc(R, b.state)
        end)
    end

    p = sortperm(xs, lt=comp)
    return (levicivita(p), RaiseLowerOps(a.ops[p]))
end

function (a::RaiseLowerOps)(X::Bases.MaybeSub{<:Bases.Slater})
    Y = deepcopy(convert(Bases.Slater, X))
    B = typeof(Y)

    sgn = 1
    for i = length(a.ops):-1:1
        sgn *= if a.ops[i] isa RaiseOp
            if a.ops[i].state in Y
                return (0, States.Zero())
            else
                create!(Y, a.ops[i].state)
            end
        else
            if !(a.ops[i].state in Y)
                return (0, States.Zero())
            else
                annihil!(Y, a.ops[i].state)
            end
        end
    end

    return (sgn, Y)
end

macro Operator(expr)
    :(@Operator($(esc(expr)), ComplexF64, DEFAULT_BASIS))
end
macro Operator(expr::Expr, basis)
    :(@Operator($(esc(expr)), ComplexF64, $(esc(basis))))
end
macro Operator(expr::Expr, T, basis)
    @assert expr.head == :(->) && expr.args[1].head == :tuple
    nbodies2 = length(expr.args[1].args)
    @assert iseven(nbodies2)
    nbodies = div(nbodies2, 2)

    quote
        x = $(esc(expr))
        FunctionOperator{$nbodies, $(esc(basis)), $(esc(T))}(x)
    end
end

Base.show(io::IO, x::RaiseOp) = print(io, "A($(x.state)')")
Base.show(io::IO, x::LowerOp) = print(io, "A($(x.state) )")
function Base.show(io::IO, x::RaiseLowerOps)
    f(a::RaiseOp) = "$(a.state)'"
    f(a::LowerOp) = "$(a.state) "

    print(io, "A(", f(x.ops[1]))
    for op in x.ops[2:end]
        print(io, ", $(f(op))")
    end
    print(io, ")")
end
@generated Base.show(io::IO, ::MIME"text/plain", x::AbstractOperator) = quote
    println(typeof(x))
    show(io, MIME"text/plain"(),
         $(if x <: FunctionOperator
               x.func
           elseif x <: ArrayOperator
               x.arr
           elseif x <: MatrixOperator
               x.mat
           else
               "Unknown AbstractOperator"
           end))
end

# Add Operator val type?
@generated Base.zero(::Type{<:Operator{N, B, <:Function, T}}) where
                    {N, B<:AbstractBasis, T} =
    :(Operator{N, B}(T, (@ntuple($(2N), _) -> zero(T))))
Base.zero(::Type{<:Operator{N, B, A}}) where {N, B<:AbstractBasis, A<:AbstractArray} =
    Operator{N, B, A}(zero(A(fill(dim(B), 2N)...)))

@generated function mask(C::Operator{N, B, A}, slots::Type{Val{S}}...) where
                        {N, S, B<:AbstractBasis, A<:AbstractArray}
    z = zero(eltype(A))
    id = one(eltype(A))

    refexpr(pos, sym) = if pos
        :(@nref($(2N), C.op, i) *= isocc(B[$sym]) ? $z : $id)
    else
        :(@nref($(2N), C.op, i) *= ~isocc(B[$sym]) ? $z : $id)
    end

    quote
        @nloops $(2N) i C.op begin
            $((refexpr(slot > 0, Symbol("i_"*string(slot))) for slot in slots)...)
        end
    end
end

# TODO
#Base.adjoint(A::Operator{N, B})

for op in (:*, :+, :-)
    @eval Base.$op(As::Operator{N, B}) where {N, B<:AbstractBasis} =
        Operator{N, B}($op(map(A -> A.op, As)...))
    @eval Base.$op(As::Operator{N, B, <:Function}) where
                  {N, B<:AbstractBasis} =
        Operator{N, B}((args...) -> $op(map(A -> A.op(args...), As)...))
end
for op in (:*, :\)
    @eval Base.$op(x::Number, A::Operator{N, B}) where {N, B<:AbstractBasis} =
        Operator{N, B}($op(x, A.op))
end
for op in (:*, :/)
    @eval Base.$op(A::Operator{N, B}, x::Number) where {N, B<:AbstractBasis} =
        Operator{N, B}($op(A.op, x))
end

end # module Operators
