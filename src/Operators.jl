@reexport module Operators
export AbstractOperator, ActionOperator, CF64ActionOperator, F64ActionOperator,
       FunctionOperator, CF64FunctionOperator, F64FunctionOperator, ArrayOperator,
       CF64ArrayOperator, F64ArrayOperator, tabulate, refop, RaiseOp, LowerOp, RaiseLowerOps,
       refop, normord

using Base.Cartesian
using Combinatorics: levicivita
using JuliaUtil: cartesian_pow
using ..Bases

abstract type AbstractOperator{B<:ConcreteBasis, T} end
struct ActionOperator{B<:ConcreteBasis, T, F<:Function} <: AbstractOperator{B, T}
    rep::F
end
const CF64ActionOperator{B<:ConcreteBasis} = ActionOperator{B, ComplexF64}
const F64ActionOperator{B<:ConcreteBasis} = ActionOperator{B, Float64}

## ?DELETEME
#ActionOperator{B, T}(f::Function) where {B<:ConcreteBasis, T} =
#    ActionOperator{B, T, typeof(f)}(f)

struct FunctionOperator{B<:ConcreteBasis, T, F<:Function} <: AbstractOperator{B, T}
    rep::F
end
const CF64FunctionOperator{B<:ConreteBasis} = FunctionOperator{B, ComplexF64}
const F64FunctionOperator{B<:ConreteBasis} = FunctionOperator{B, Float64}

## ?DELETEME
#FunctionOperator{B, T}(f::Function) where {B<:ConcreteBasis, T} =
#    FunctionOperator{B, T, typeof(f)}(f)
#
struct ArrayOperator{B<:ConcreteBasis, T, A<:AbstractArray{T}} <: AbstractOperator{B, T}
    rep::A

    function ArrayOperator{B, T, A}(rep::A) where {B<:ConcreteBasis, T, A<:AbstractArray{T}}
        @assert size(rep) == (innerdims(B)..., innerdims(B))
        new(rep)
    end
end
const CF64ArrayOperator{B<:ConcreteBasis} = ArrayOperator{B, ComplexF64, <:Array{ComplexF64}}
const F64ArrayOperator{B<:ConcreteBasis} = ArrayOperator{N, B, Float64, <:Array{Float64}}

## ?DELETEME
#ArrayOperator{B}(a::AbstractArray) where B<:ConcreteBasis =
#    ArrayOperator{B, eltype(a), typeof(a)}(a)

#=ManyBody.=#rep(op::Union{ActionOperator, FunctionOperator, ArrayOperator}) = op.rep
#ManyBody.reptype(op::AbstractOperator) = reptype(typeof(op))
#=ManyBody.=#reptype(::Type{<:ArrayOperator{<:Any, <:Any, A}}) where A = A
#=ManyBody.=#reptype(::Type{<:ActionOperator{<:Any, <:Any, F}}) where F = F
#=ManyBody.=#reptype(::Type{<:FunctionOperator{<:Any, <:Any, F}}) where F = F
Base.eltype(::Type{<:AbstractOperator{<:Any, T}}) where T = T

#ManyBody.basistype(x::AbstractOperator) = basistype(typeof(x))
ManyBody.basistype(::Type{<:AbstractOperator{B}}) where B<:ConcreteBasis = B
Bases.rank(O::Type{<:AbstractOperator}) = rank(basistype(O))

### Definition
## For now, this does not work
#(op::AbstractOperator{N, <:Any, T})(args...) where {N, T} = op(Array{T, N}, args...)
## Work around
(op::ActionOperator)(args...) = applyop(op, args...)
(op::FunctionOperator)(args...) = applyop(op, args...)
(op::ArrayOperator)(args...) = applyop(op, args...)

### Kernels
## Subsumes FunctionOperator{} method
applyop(op::AbstractOperator{B}, arg::B) where B<:ConcreteBasis = sum(B) do b
    op[b, arg]*b
end
applyop(op::ActionOperator{B}, arg::B) where B<:ConcreteBasis = rep(op)(arg)
applyop(op::ArrayOperator{B}, arg::Bases.MaybeIndex{B}) where B<:ConcreteBasis =
    sum(eachindex(B)) do I
        op[b, index(arg)]*convert(B, arg)
    end

### Dispatch
applyop(op::AbstractOperator, arg) = applyop(op, convert(basistype(op), arg))
applyop(op::AbstractOperator, arg::Bases.Neg) = -applyop(op, inner(arg))
applyop(op::AbstractOperator, arg::States.Scaled) = arg.coeff*applyop(op, arg.state)
function applyop(op::AbstractOperator, arg::ArrayState)
    OB, AB = basistype.((op, arg))
    OBI, ABI = indextype.((OB, AB))
    T = promote_type(eltype(op), eltype(arg))
    ret = similar(reptype(arg), T, innerdims(OB)) |> zero

    # Assumes findall() always returns the correct index type...
    for I in findall(!iszero, rep(arg))
        ret[OBI[ABI[I]]] = rep(arg)[I]*applyop(op, AB[I])
    end
end

### Definition
Base.getindex(op::AbstractOperator, args...) = matrixelem(op, args...)

### Kernels
## Subsumes ActionOperator{} method
matrixelem(op::AbstractOperator{B}, argl::B, argr::B) where B<:ConcreteBasis =
    convert(eltype(op), argl'op(argr))
matrixelem(op::FunctionOperator{B}, argl::B, argr::B) where B<:ConcreteBasis =
    convert(eltype(op), rep(op)(argl, argr))
matrixelem(op::ArrayOperator{B}, argl::Bases.MaybeIndex{B}, argr::Bases.MaybeIndex{B}) where
          B<:ConcreteBasis =
    rep(op)[index(argl), index(argr)]

### Dispatch
function matrixelem(op::AbstractOperator, args...)  
    N = rank(op)
    matrixelem(op, prod(args[1:N]), prod(args[N+1:end]))
end
function matrixelem(op::AbstractOperator, argl, argr)
    B = basistype(op)
    matrixelem(op, convert(B, argl), convert(B, argr))
end
matrixelem(op::AbstractOperator, argl::Bases.Neg, argr::Bases.Neg) =
    matrixeleme(op, inner(argl), inner(argr))
@commutes (2,3) conj matrixelem(op::AbstractOperator, argl::Bases.Neg, argr) =
    -matrixelem(op, inner(argl), argr)
matrixelem(op::AbstractOperator, argl::States.Scaled, argr::States.Scaled) =
    conj(argl.coeff)*argr.coeff*matrixelem(op, argl.state, argr.state)
@commutes (2,3) conj matrixelem(op::AbstractOperator, argl, argr::States.Scaled) =
    argr.coeff*matrixelem(op, argl, argr.state)
## Add special case for ArrayOperator{}?
function matrixelem(op::AbstractOperator, argl::ArrayState, argr::ArrayState)
    Bl, Br = basistype.((argl, argr))

    sum(Iterators.product(eachindex(Bl), eachindex(Br))) do K
        I, J = K
        
        l, r = (Bl[I], Br[J])
        matrixelem(op, l, r)*(argl'l)*(r'argr)
    end
end

### Update Line
##############################################################################################

nbodies(::Type{<:AbstractOperator{N}}) where N = N
nbodies(op::AbstractOperator) = nbodies(typeof(op))

tabulate(op) = tabulate(Array{ComplexF64}, op)
tabulate(::Type{A}, x::Number) where A<:AbstractArray = convert(eltype(A), x)
tabulate(::Type{A}, x::AbstractArray) where A<:AbstractArray = convert(A, x)
tabulate(::Type{A}, op::ArrayOperator{N, B}) where {A<:AbstractArray, N, B<:Bases.Index} =
    ArrayOperator{N, B}(convert(A, rep(op)))
tabulate(::Type{A}, op::ArrayOperator{N, B}) where {A<:AbstractArray, N, B<:AbstractBasis} =
    ArrayOperator{N, Bases.Index{B}}(convert(A, rep(op)))
#@generated function tabulate(::Type{A}, op::AbstractOperator{N, B}) where
#                            {A<:AbstractArray, N, B<:AbstractBasis}
#    quote
#        arr = similar(A, $(fill(dim(B), 2N)...))
#        @nloops $(2N) s (_ -> B) begin
#            @nref($(2N), arr, i -> index(s_i)) = @nref($(2N), op, s)
#        end
#
#        ArrayOperator{N, Bases.Index{B}}(arr)
#    end
#end
function tabulate(::Type{A}, op::AbstractOperator) where {A<:AbstractArray}
    arr = similar(A, fill(dim(basistype(op)), 2nbodies(op))...)
    for S in cartesian_pow(basistype(op), Val{2nbodies(op)})
        arr[CartesianIndex(map(index, S))] = op[S...]
    end

    ArrayOperator{nbodies(op), Bases.Index{basistype(op)}}(arr)
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

    (levicivita(p), RaiseLowerOps(a.ops[p]))
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
@generated Base.zero(::Type{<:ActionOperator{N, B, T}}) where
                    {N, B<:AbstractBasis, T} =
    :(ActionOperator{N, B, T}((@ntuple($N, _) -> zeros(T, $(fill(dim(B), N)...)))))
@generated Base.zero(::Type{<:FunctionOperator{N, B, T}}) where
                    {N, B<:AbstractBasis, T} =
    :(FunctionOperator{N, B, T}((@ntuple($(2N), _) -> zero(T))))
@generated Base.zero(::Type{<:ArrayOperator{N, B, <:Any, A}}) where
                    {N, B<:AbstractBasis, A<:AbstractArray} =
    :(ArrayOperator{N, B}(zero(similar(A, $(fill(dim(B), 2N)...)))))

#@generated Base.adjoint(op::ActionOperator{N, B, T}) where {N, B<:AbstractBasis, T} =
#    :(ActionOperator{N, B, T}(@ntuple($N, b) ->
#                                    conj.(permutedims(@ncall($N, rep(op), b), $N:-1:1))))
@generated function Base.adjoint(op::FunctionOperator{N, B, T}) where {N, B<:AbstractBasis, T}
    reversed = [Symbol("arg_$i") for i = 2N:-1:1]
    :(FunctionOperator{N, B, T}(@ntuple($(2N), arg) -> conj(rep(op)($(reversed...)))))
end
Base.adjoint(op::ArrayOperator{N, B}) where {N, B<:AbstractBasis} =
    ArrayOperator{N, B}(conj.(permutedims(rep(op), 2N:-1:1)))

for op in (:*, :+, :-)
    @eval Base.$op(As::ArrayOperator{N, B}...) where {N, B<:AbstractBasis} =
        ArrayOperator{N, B}($op(map(A -> rep(A), As)...))
    ## Should allow different T's. Also, @generate this
    @eval Base.$op(As::ActionOperator{N, B, T}...) where {N, B<:AbstractBasis, T} =
        ActionOperator{N, B, T}((args...) -> $op(map(A -> rep(A)(args...), As)...))
    @eval Base.$op(As::FunctionOperator{N, B, T}...) where {N, B<:AbstractBasis, T} =
        FunctionOperator{N, B, T}((args...) -> $op(map(A -> rep(A)(args...), As)...))
end
for op in (:*, :\)
    @eval Base.$op(x::Number, A::ArrayOperator{N, B}) where {N, B<:AbstractBasis} =
        ArrayOperator{N, B}($op(x, rep(A)))
    ## Should allow different T's...
    @eval Base.$op(x::Number, A::ActionOperator{N, B, T}) where {N, B<:AbstractBasis, T} =
        ActionOperator{N, B, T}((args...) -> $op(x, rep(op)(args...)))
    @eval Base.$op(x::Number, A::FunctionOperator{N, B, T}) where {N, B<:AbstractBasis, T} =
        FunctionOperator{N, B, T}((args...) -> $op(x, rep(op)(args...)))
end
for op in (:*, :/)
    @eval Base.$op(A::ArrayOperator{N, B}, x::Number) where {N, B<:AbstractBasis} =
        ArrayOperator{N, B}($op(rep(A), x))
    ## Should allow different T's...
    @eval Base.$op(A::ActionOperator{N, B, T}, x::Number) where {N, B<:AbstractBasis, T} =
        ActionOperator{N, B, T}((args...) -> $op(rep(A)(args...), x))
    @eval Base.$op(A::FunctionOperator{N, B, T}, x::Number) where {N, B<:AbstractBasis, T} =
        FunctionOperator{N, B, T}((args...) -> $op(rep(A)(args...), x))
end

end # module Operators
