@reexport module Operators
export Operator, tabulate, refop, RaiseOp, LowerOp, RaiseLowerOps, A, refop, contract,
       normord, apply_normord_rl, @Operator

using Combinatorics: levicivita
using Loppy.Util: cartesian_pow
using ..Bases
using ..States

struct Operator{N, B<:AbstractBasis, Op, T}
    op::Op

    Operator{N, B, Op, T}(op::Op) where {N, B, Op, T} = new(op)
    function Operator{N, B, Op, T}(op::Op) where {N, B, T, Op<:AbstractArray}
        @assert eltype(Op) == T
        new(op)
    end
end
const COperator{N, B<:AbstractBasis, Op} = Operator{N, B, Op, Complex64}
Operator{N, B}(T, op) where {N, B<:AbstractBasis} = Operator{N, B, typeof(op), T}(op)
Operator{N, B}(op::AbstractArray{T}) where {N, B<:AbstractBasis} =
    Operator{N, B, typeof(op), T}(op)

@generated (op::Operator{N, B, Op, T})(args::Vararg{<:Bases.MaybeSub{B}, N2}) where
                                   {N, N2, B<:AbstractBasis, Op, T} =
    :(applyop(Val{$(2N)}, op, $((:(args[$i]) for i=1:N2)...)))
@generated applyop(::Type{Val{N2}}, op::Operator{N, B, <:Function, T},
                   sps::Vararg{<:Bases.MaybeSub{B}, N2}) where {N, N2, B<:AbstractBasis, T} =
    :(convert(T, op.op($((sps[i] !== B && sps[i] <: Bases.Sub? :(sps[$i].state) : :(sps[$i])
                          for i = 1:N2)...))))
@generated applyop(::Type{Val{N2}}, op::Operator{N, B, <:AbstractArray{T, N2}},
                   sps::Vararg{Union{B, SB}, N2}) where
                  {N, N2, B<:AbstractBasis, T, SB<:Bases.Sub{B}} =
    :(op.op[$((sps[i] !== B && <: Bases.Sub? :(index(sps[$i].state)) : :(index(sps[$i]))
               for i = 1:N2)...)])

nbodies(::Type{<:Operator{N}}) where N = N
nbodies(op::Operator) = nbodies(typeof(op))

tabulate(op) = tabulate(Complex64, op)
tabulate(::Type, x::Number) = x
tabulate(::Type, x::AbstractArray) = x
tabulate(::Type{T}, op::Operator{N, B, <:AbstractArray}) where {T, N, B<:Bases.Index} =
    Operator{N, B}(convert(AbstractArray{T}, op.op))
function tabulate(::Type{T}, op::Operator{N, B}) where {T, N, B<:AbstractBasis}
    mat = Array{T}(fill(dim(B), 2N)...)
    for ss in cartesian_pow(B, 2N)
        mat[CartesianIndex(map(index, ss))] = op(ss...)
    end

    Operator{N, Bases.Index{B}}(mat)
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
Base.ctranspose(op::RaiseOp{B}) where B<:AbstractBasis = LowerOp{B}(op.state)
Base.ctranspose(op::LowerOp{B}) where B<:AbstractBasis = RaiseOp{B}(op.state)
Base.ctranspose(op::RaiseLowerOps{B}) where B<:AbstractBasis =
    RaiseLowerOps{B}(reverse(map(ctranspose, op.ops)))

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
    :(@Operator($(esc(expr)), Complex64, DEFAULT_BASIS))
end
macro Operator(expr::Expr, basis)
    :(@Operator($(esc(expr)), Complex64, $(esc(basis))))
end
macro Operator(expr::Expr, T, basis)
    @assert expr.head == :(->) && expr.args[1].head == :tuple
    nbodies2 = length(expr.args[1].args)
    @assert iseven(nbodies2)
    nbodies = div(nbodies2, 2)

    quote
        x = $(esc(expr))
        Operator{$nbodies, $(esc(basis))}($(esc(T)), x)
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

# Add Operator val type?
@generated Base.zero(::Type{<:Operator{N, B, <:Function, T}}) where
                    {N, B<:AbstractBasis, T} =
    :(Operator{N, B}(T, (@ntuple($(2N), _) -> zero(T)))
Base.zero(::Type{<:Operator{N, B, A}}) where A<:AbstractArray =
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
Base.ctranspose(A::Operator{N, B})

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
