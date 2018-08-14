@reexport module Operators
export AbstractOperator, FunctionOperator, CFunctionOperator, ArrayOperator, tabulate, refop,
       RaiseOp, LowerOp, RaiseLowerOps, A, refop, contract, normord 

using Base.Cartesian
using Combinatorics: levicivita
using JuliaUtil: cartesian_pow
using ..Bases

abstract type AbstractOperator{N, B<:AbstractBasis, T} end
struct FunctionOperator{N, B<:AbstractBasis, T, F<:Function} #=
    =# <: AbstractOperator{N, B, T}
    rep::F
end
const CF64FunctionOperator{N, B<:AbstractBasis} = FunctionOperator{N, B, ComplexF64}
const F64FunctionOperator{N, B<:AbstractBasis} = FunctionOperator{N, B, Float64}

FunctionOperator{N, B, T}(f::Function) where {N, B<:AbstractBasis, T} =
    FunctionOperator{N, B, T, typeof(f)}(f)

struct ArrayOperator{N, B<:AbstractBasis, T, A<:AbstractArray{T}} #=
    =# <: AbstractOperator{N, B, T}
    rep::A

    function ArrayOperator{N, B, T, A}(rep::A) where
                                      {N, B<:AbstractBasis, T, A<:AbstractArray{T}}
        @assert ndims(rep) == 2N
        new(rep)
    end
end
const CF64ArrayOperator{N, B<:AbstractBasis} =
    ArrayOperator{N, B, ComplexF64, <:Array{ComplexF64}}
const F64ArrayOperator{N, B<:AbstractBasis} =
    ArrayOperator{N, B, Float64, <:Array{Float64}}

function ArrayOperator{N, B, T, A}(a::AbstractArray) where
                                  {N, B<:AbstractBasis, T, A<:AbstractArray{T}}
    a = convert(A, a)
    ArrayOperator{N, B, T, A}(a)
end
function ArrayOperator{N, B, T}(a::AbstractArray) where {N, B<:AbstractBasis, T}
    a = convert(AbstractArray{T}, a)
    ArrayOperator{N, B, T, typeof(a)}(a)
end
ArrayOperator{N, B}(a::AbstractArray) where {N, B<:AbstractBasis} =
    ArrayOperator{N, B, eltype(a), typeof(a)}(a)
function ArrayOperator(B::AbstractBasis, a::AbstractArray)
    @assert iseven(ndims(a))
    ArrayOperator{div(ndims(a), 2), B, eltype(a), typeof(a)}(a)
end

rep(op::AbstractOperator) = op.rep
reptype(op::ArrayOperator{<:Any, <:Any, <:Any, A}) where A = A
reptype(op::FunctionOperator{<:Any, <:Any, <:Any, F}) where F = F
Base.eltype(::Type{<:AbstractOperator{<:Any, <:Any, T}}) where T = T

## For now, this does not work
#(op::AbstractOperator{N, <:Any, T})(args...) where {N, T} = op(Array{T, N}, args...)
## Work around
(op::FunctionOperator{N, <:Any, T})(args...) where {N, T} = op(Array{T, N}, args...)
(op::ArrayOperator{N, <:Any, T})(args...) where {N, T} = op(Array{T, N}, args...)
@generated (op::FunctionOperator{N, B, T})(args::Vararg{B, N}) where
                                          {N, B<:AbstractBasis, T, A<:AbstractArray{T, N}} =
    :(@ncall($N, rep(op), i -> args[i]))
(op::FunctionOperator{N, B, T})(::Type{A}, args::Vararg{B, N}) where
                               {N, B<:AbstractBasis, T, A<:AbstractArray{T, N}} =
    convert(A, op(args...))
@generated (op::FunctionOperator{N, B})(::Type{A}, state::AbstractArray{<:Any, N}) where
                                       {N, B<:AbstractBasis, T, A<:AbstractArray{T, N}} = quote
    ret = similar(A, $(fill(dim(B), N)...))
    @nloops $N b (_ -> B) begin
        @nref($N, ret, i -> index(b_i)) = @nref($N, state, i -> index(b_i))*@ncall($N, rep(op), b)
    end

    ret
end

@generated function Base.getindex(op::AbstractOperator{N, B}, args::Vararg{B, N2}) where
                                 {N, N2, B<:AbstractBasis}
    @assert N2 == 2N
    :(@ncall($N2, matrixelem, op, i -> args[i]))
end

@generated matrixelem(op::ArrayOperator{N, B}, args::Vararg{B, N2}) where
                     {N, N2, B<:AbstractBasis} =
    :(@nref($N2, rep(op), i -> index(args[i])))
@generated matrixelem(op::AbstractOperator{N, B, T}, args::Vararg{B, N2}) where
                     {N, N2, T, B<:AbstractBasis} =
    :(convert(T, prod(args[1:$N])'*@ncall($N, op, i -> args[$N+i])))

nbodies(::Type{<:AbstractOperator{N}}) where N = N
nbodies(op::AbstractOperator) = nbodies(typeof(op))

tabulate(op) = tabulate(Array{ComplexF64}, op)
tabulate(::Type{A}, x::Number) where A<:AbstractArray = convert(eltype(A), x)
tabulate(::Type{A}, x::AbstractArray) where A<:AbstractArray = convert(A, x)
tabulate(::Type{A}, op::ArrayOperator{N, B}) where {A<:AbstractArray, N, B<:Bases.Index} =
    ArrayOperator{N, B}(convert(A, rep(op)))
tabulate(::Type{A}, op::ArrayOperator{N, B}) where {A<:AbstractArray, N, B<:AbstractBasis} =
    ArrayOperator{N, Bases.Index{B}}(convert(A, rep(op)))
@generated function tabulate(::Type{A}, op::AbstractOperator{N, B}) where
                            {A<:AbstractArray, N, B<:AbstractBasis}
    quote
        arr = similar(A, $(fill(dim(B), 2N)...))
        @nloops $(2N) s (_ -> B) begin
            @nref($(2N), arr, i -> index(s_i)) = @nref($(2N), op, s)
        end

        ArrayOperator{N, Bases.Index{B}}(arr)
    end
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

@generated (a::RaiseLowerOps)(X::Bases.MaybeSub{<:Bases.Slater}) = quote
    z = zero(Vector(X))
    Y = deepcopy(convert(Bases.Slater, X))
    B = typeof(Y)

    sgn = 1
    for i = length(a.ops):-1:1
        sgn *= if a.ops[i] isa RaiseOp
            if a.ops[i].state in Y
                return ZeroState()
            else
                create!(Y, a.ops[i].state)
            end
        else
            if !(a.ops[i].state in Y)
                return ZeroState()
            else
                annihil!(Y, a.ops[i].state)
            end
        end
    end

    # Probably a better way to do this.
    $(if X <: Bases.Sub
          quote
              Z = convert(typeof(X), Y)
              if Z in typeof(X)
                  Y = Z
              end
          end
      else
          :()
      end)
    sgn == 1 ? Y : sgn*Y
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
function Base.show(io::IO, mime::MIME"text/plain", x::AbstractOperator)
    println(typeof(x))
    show(io, mime, rep(x))
end

# Add Operator val type?
@generated Base.zero(::Type{<:FunctionOperator{N, B, T}}) where
                    {N, B<:AbstractBasis, T} =
    :(FunctionOperator{N, B, T}((@ntuple($N, _) -> zeros(T, $(fill(dim(B), N)...)))))
@generated Base.zero(::Type{<:ArrayOperator{N, B, <:Any, A}}) where
                    {N, B<:AbstractBasis, A<:AbstractArray} =
    :(ArrayOperator{N, B}(zero(similar(A, $(fill(dim(B), 2N)...)))))

#@generated Base.adjoint(op::FunctionOperator{N, B, T}) where {N, B<:AbstractBasis, T} =
#    :(FunctionOperator{N, B, T}(@ntuple($N, b) ->
#                                    conj.(permutedims(@ncall($N, rep(op), b), $N:-1:1))))
Base.adjoint(op::ArrayOperator{N, B}) where {N, B<:AbstractBasis} =
    ArrayOperator{N, B}(conj.(permutedims(rep(op), 2N:-1:1)))

for op in (:*, :+, :-)
    @eval Base.$op(As::ArrayOperator{N, B}...) where {N, B<:AbstractBasis} =
        ArrayOperator{N, B}($op(map(A -> rep(A), As)...))
    ## Should allow different T's...
    @eval Base.$op(As::FunctionOperator{N, B, T}...) where {N, B<:AbstractBasis, T} =
        FunctionOperator{N, B, T}((args...) -> $op(map(A -> rep(A)(args...), As)...))
end
for op in (:*, :\)
    @eval Base.$op(x::Number, A::ArrayOperator{N, B}) where {N, B<:AbstractBasis} =
        ArrayOperator{N, B}($op(x, rep(A)))
    ## Should allow different T's...
    @eval Base.$op(x::Number, A::FunctionOperator{N, B, T}) where {N, B<:AbstractBasis, T} =
        FunctionOperator{N, B, T}((args...) -> $op(x, rep(op)(args...)))
        
end
for op in (:*, :/)
    @eval Base.$op(A::ArrayOperator{N, B}, x::Number) where {N, B<:AbstractBasis} =
        ArrayOperator{N, B}($op(rep(A), x))
    ## Should allow different T's...
    @eval Base.$op(A::FunctionOperator{N, B, T}, x::Number) where {N, B<:AbstractBasis, T} =
        FunctionOperator{N, B, T}((args...) -> $op(rep(A)(args...), x))
end

end # module Operators
