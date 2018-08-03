@reexport module Operators
export Operator, tabulate, refop, RaiseOp, LowerOp, RaiseLowerOps, A, refop, contract,
       normord, apply_normord_rl, @Operator#, @OperatorNBody

using Combinatorics: permutations, levicivita
using Loppy.Util: cartesian_pow
using ..Bases
using ..States

struct Operator{N, B<:AbstractBasis, Op}
    op::Op
end
@generated (op::Operator{N, B, Op})(args::Vararg{<:Bases.MaybeSub{B}, N2}) where {N, N2, B, Op} =
    :(applyop(Val{$(2N)}, op, $((:(args[$i]) for i=1:N2)...)))
@generated applyop(::Type{Val{N2}}, op::Operator{N, B, <:Function},
                   sps::Vararg{<:Bases.MaybeSub{B}, N2}) where {N, N2, B} =
    :(op.op($((sps[i] !== B && sps[i] <: Bases.Sub? :(sps[$i].state) : :(sps[$i])
               for i = 1:N2)...)))
@generated applyop(::Type{Val{N2}}, op::Operator{N, B, <:AbstractArray{T, N2}},
                   sps::Vararg{Union{B, SB}, N2}) where {N, N2, B, T, SB<:Bases.Sub{B}} =
    :(op.op[$((sps[i] !== B && <: Bases.Sub? :(index(sps[$i].state)) : :(index(sps[$i]))
               for i = 1:N2)...)])

nbodies(::Type{<:Operator{N}}) where N = N
nbodies(op::Operator) = nbodies(typeof(op))

tabulate(op) = tabulate(Complex64, op)
tabulate(op::Operator{N, <:Bases.Index, <:AbstractArray}) where N = op
function tabulate(::Type{T}, op::Operator{N, B, Op}) where {T, N, B, Op}
    mat = Array{T}(fill(dim(B), 2N)...)
    for ss in cartesian_pow(B, 2N)
        mat[CartesianIndex(map(index, ss))] = op(ss...)
    end

    Operator{N, Bases.Index{B}, typeof(mat)}(mat)
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

Base.convert(::Type{RaiseLowerOps{B}}, a::RLOp{B}) where B = RaiseLowerOps{B}([a])
Base.convert(::Type{RaiseLowerOps}, a::RLOp{B}) where B = RaiseLowerOps{B}([a])

Base.:*(op1::RLOp{B}, op2::RLOp{B}) where B = RaiseLowerOps{B}([op1, op2])
Base.:*(op1::RLOp{B}, op2::RaiseLowerOps{B}) where B = RaiseLowerOps{B}([op1; op2.ops])
Base.:*(op1::RaiseLowerOps{B}, op2::RLOp{B}) where B = RaiseLowerOps{B}([op2.ops; op2])
Base.:*(op1::RaiseLowerOps{B}, op2::RaiseLowerOps{B}) where B =
    RaiseLowerOps{B}([op1.ops; op2.ops])
Base.ctranspose(op::RaiseOp{B}) where B = LowerOp{B}(op.state)
Base.ctranspose(op::LowerOp{B}) where B = RaiseOp{B}(op.state)
Base.ctranspose(op::RaiseLowerOps{B}) where B =
    RaiseLowerOps{B}(reverse(map(ctranspose, op.ops)))

A(sps...) = A(sps)
A(T::Type, sps...) = A(T, sps)

A(p::Bases.MaybeSub) = LowerOp(p)
A(p::Bra{<:Bases.MaybeSub}) = RaiseOp(p.state)

@generated function A(sps::NTuple{N, Union{T, Bra{T}, Raised{T}}}) where
                     {B, T<:Bases.MaybeSub{B}, N}
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

    :(RaiseLowerOps{B}([$(args...)]))
end

n_ops(::RLOp) = 1
n_ops(a::RaiseLowerOps) = length(a.ops)

refop(s::Bases.Sub{B}) where B<:AbstractBasis = refop(s.state)
function refop(s::Bases.Slater{B}) where B
    a = RaiseLowerOps{B}(map(x -> RaiseOp{B}(indexp(R, x)), find(s.parts)))
    b = RaiseLowerOps{B}(map(x -> LowerOp{B}(indexh(R, x)), find(s.holes)))
    a*b
end

contract(::Type, a1::OP, a2::OP) where OP<:RLOp = 0
contract(::Type{R}, a1::RaiseOp{B}, a2::LowerOp{B}) where {B, R<:RefState{B}} =
    a1.state == a2.state ? isocc(R, a1.state) : 0
contract(::Type{R}, a1::LowerOp{B}, a2::RaiseOp{B}) where {B, R<:RefState{B}} =
    a1.state == a2.state ? ~isocc(R, a1.state) : 0
contract(::Type, a::RLOp) = 0

function contract(::Type{R}, a::RaiseLowerOps{B}) where {B, R<:RefState{B}}
    numops = n_ops(a)
    isodd(numops) && return 0

    sum(permutations(1:numops)) do perm
        prod(1:div(numops, 2)) do k
            levicivita(perm)*contract(B, a.ops[perm[2k-1]], a.ops[perm[2k]])
        end
    end
end

normord(a::RaiseLowerOps{B}) where B = normord(RefStates.Vacuum{B}, a)
function normord(::Type{R}, a::RaiseLowerOps{B}) where {B, R<:RefState{B}}
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

function apply_normord_rl(a::RaiseLowerOps, X::Bases.MaybeSub{<:Bases.Slater})
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

end # module Operators
