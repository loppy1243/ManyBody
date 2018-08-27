@reexport module Operators
export AbstractOperator, ActionOperator, CF64ActionOperator, F64ActionOperator,
       FunctionOperator, CF64FunctionOperator, F64FunctionOperator, ArrayOperator,
       CF64ArrayOperator, F64ArrayOperator, tabulate, refop, RaiseOp, LowerOp, RaiseLowerOp,
       refop, normord, matrixiter
using ..ManyBody
using Base.Cartesian
using ..@commutes, ..@disallow
using Combinatorics: levicivita
using JuliaUtil: cartesian_pow

abstract type AbstractOperator{B<:ConcreteBasis, T} end
struct ActionOperator{B<:ConcreteBasis, T, F<:Function} <: AbstractOperator{B, T}
    rep::F
end
const CF64ActionOperator{B<:ConcreteBasis} = ActionOperator{B, ComplexF64}
const F64ActionOperator{B<:ConcreteBasis} = ActionOperator{B, Float64}

ActionOperator{B, T}(f::Function) where {B<:ConcreteBasis, T} =
    ActionOperator{B, T, typeof(f)}(f)

## ?DELETEME
#ActionOperator{B, T}(f::Function) where {B<:ConcreteBasis, T} =
#    ActionOperator{B, T, typeof(f)}(f)

struct FunctionOperator{B<:ConcreteBasis, T, F<:Function} <: AbstractOperator{B, T}
    rep::F
end
const CF64FunctionOperator{B<:ConcreteBasis} = FunctionOperator{B, ComplexF64}
const F64FunctionOperator{B<:ConcreteBasis} = FunctionOperator{B, Float64}

FunctionOperator{B, T}(f::Function) where {B<:ConcreteBasis, T} =
    FunctionOperator{B, T, typeof(f)}(f)

## ?DELETEME
#FunctionOperator{B, T}(f::Function) where {B<:ConcreteBasis, T} =
#    FunctionOperator{B, T, typeof(f)}(f)
#
struct ArrayOperator{B<:ConcreteBasis, T, A<:AbstractArray{T}} <: AbstractOperator{B, T}
    rep::A

    function ArrayOperator{B, T, A}(rep::A) where {B<:ConcreteBasis, T, A<:AbstractArray{T}}
        @assert size(rep) == (innerdims(B)..., innerdims(B)...)
        new(rep)
    end
end
const CF64ArrayOperator{B<:ConcreteBasis, N2} =  
    ArrayOperator{B, ComplexF64, Array{ComplexF64, N2}}
const F64ArrayOperator{B<:ConcreteBasis, N2} = ArrayOperator{B, Float64, Array{Float64, N2}}

ArrayOperator{B}(rep::AbstractArray) where B<:ConcreteBasis =
    ArrayOperator{B, eltype(rep), typeof(rep)}(rep)

## ?DELETEME
#ArrayOperator{B}(a::AbstractArray) where B<:ConcreteBasis =
#    ArrayOperator{B, eltype(a), typeof(a)}(a)

ManyBody.reptype(::Type{<:ArrayOperator{<:Any, <:Any, A}}) where A = A
ManyBody.reptype(::Type{<:ActionOperator{<:Any, <:Any, F}}) where F = F
ManyBody.reptype(::Type{<:FunctionOperator{<:Any, <:Any, F}}) where F = F
Base.eltype(::Type{<:AbstractOperator{<:Any, T}}) where T = T

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
applyop!(op::AbstractOperator, ::States.Zero) = States.Zero()
applyop(op::AbstractOperator, ::States.Zero) = States.Zero()
## Subsumes FunctionOperator{} method
applyop(op::AbstractOperator{B}, arg::B) where B<:ConcreteBasis = sum(B) do b
    op[b, arg]*b
end
applyop(op::ActionOperator{B}, arg::B) where B<:ConcreteBasis = op.rep(arg)
applyop(op::ArrayOperator{B}, arg::Bases.MaybeIndex{B}) where B<:ConcreteBasis =
    sum(eachindex(B)) do I
        op[b, index(arg)]*convert(B, arg)
    end

### Dispatch
applyop(op::AbstractOperator, arg) = applyop(op, convert(basistype(op), arg))
applyop!(op::AbstractOperator, arg::Bases.Neg) = -applyop!(op, inner(arg))
applyop(op::AbstractOperator, arg::Bases.Neg) = -applyop(op, inner(arg))
applyop!(op::AbstractOperator, arg::States.Scaled) = arg.coeff*applyop!(op, arg.state)
applyop(op::AbstractOperator, arg::States.Scaled) = arg.coeff*applyop(op, arg.state)
function applyop!(op::AbstractOperator{B}, arg::ArrayState{B}) where B<:ConcreteBasis
    ret = applyop(op, arg)
    rep(arg) .= rep(convert(ArrayState, ret))

    arg
end
function applyop(op::AbstractOperator, arg::ArrayState)
    OB, AB = basistype.((op, arg))
    OBI, ABI = indextype.((OB, AB))
    T = promote_type(eltype(op), eltype(arg))
    ret = similar(reptype(arg), T, innerdims(OB)) |> zero

    # Assumes findall() always returns the correct index type...
    for I in findall(!iszero, arg.rep)
        ret += arg.rep[I]*applyop(op, AB[I])
    end

    ret
end

### Definition
Base.getindex(op::AbstractOperator, args...) = matrixelem(op, args...)

### Kernels
## Subsumes ActionOperator{} method
matrixelem(op::AbstractOperator{B}, argl::B, argr::B) where B<:ConcreteBasis =
    convert(eltype(op), argl'op(argr))

matrixelem(op::FunctionOperator{B}, argl::B, argr::B) where B<:ConcreteBasis =
    convert(eltype(op), op.rep(argl, argr))

matrixelem(op::ArrayOperator{B}, argl::B, argr::B) where B<:ConcreteBasis =
    op.rep[index(argl), index(argr)]
matrixelem(op::ArrayOperator{B}, argl::Bases.MaybeIndex{B}, argr::Bases.MaybeIndex{B}) where
          B<:ConcreteBasis =
    op.rep[index(argl), index(argr)]

### Dispatch
function matrixelem(op::AbstractOperator, args...)  
    N = rank(op)
    matrixelem(op, prod(args[1:N]), prod(args[N+1:end]))
end
@generated function matrixelem(op::AbstractOperator, args::Union{Int, Base.CartesianIndex}...)
    N = rank(op)
    B = basistype(op)
    BI = indextype(B)

    exprs = map(enumerate(args)) do X
        i, T = X
        if T <: Int
            :(args[$i])
        else
            Expr(:..., :(Tuple(args[$i])))
        end
    end
    tuple_expr = :(($(exprs...),))
    quote
        x = $tuple_expr
        matrixelem(op, $BI(x[1:$N]), $BI(x[$N+1:end]))
    end
end
function matrixelem(op::AbstractOperator, argl::AbstractBasis, argr::AbstractBasis)
    B = basistype(op)
    matrixelem(op, convert(B, argl), convert(B, argr))
end

### I don't think we want these, as index notation should be reserved for the actual matrix
### elements, and off-basis matrix elements should go through the action notation.

#matrixelem(op::AbstractOperator, argl, argr) = matrixelem(op, promote(argl, argr)...)

#matrixelem(op::AbstractOperator, argl::Bases.Neg, argr::Bases.Neg) =
#    matrixelem(op, inner(argl), inner(argr))
#@commutes (2,3) conj matrixelem(op::AbstractOperator, argl::Bases.Neg, argr) =
#    -matrixelem(op, inner(argl), argr)
#
#matrixelem(op::AbstractOperator, argl::States.Scaled, argr::States.Scaled) =
#    conj(argl.coeff)*argr.coeff*matrixelem(op, argl.state, argr.state)
#@commutes (2,3) conj matrixelem(op::AbstractOperator, argl, argr::States.Scaled) =
#    argr.coeff*matrixelem(op, argl, argr.state)
#
#function matrixelem(op::AbstractOperator, argl::ArrayState, argr::ArrayState)
#    BlI, BrI = indextype.(basistype.((argl, argr)))
#
#    sum(Iterators.product(eachindex(Bl), eachindex(Br))) do K
#        I, J = K
#        
#        L, R = (BlI[I], BrI[J])
#        matrixelem(op, L, R)*(argl'L)*(R'argr)
#    end
#end

matrixiter(op::AbstractOperator) = cartesian_pow(basistype(op), Val{2})

tabulate(op) = tabulate(Array{eltype(op)}, op)
tabulate(A::Type{<:AbstractArray}, op::ArrayOperator) =
    ArrayOperator{basistype(op)}(convert(A, op.rep))
function tabulate(A::Type{<:AbstractArray}, op::AbstractOperator)
    B = basistype(op)
    BI = indextype(B)
    ds = innerdims(B)
    ret = similar(A, ds..., ds...)

    for (I, J) in cartesian_pow(eachindex(B), Val{2})
        ret[I, J] = op[BI[I], BI[J]]
    end

    ArrayOperator{B}(ret)
end


### Could make RaiseLowerOp{} type stable by combining RaiseOp{} and LowerOp{} into one type
### with a parameter to select whether it is a raise or lower op.

struct RaiseOp{SP<:ConcreteBasis} <: AbstractOperator{Bases.Slater{SP}, Int}
    state::SP
end
struct LowerOp{SP<:ConcreteBasis} <: AbstractOperator{Bases.Slater{SP}, Int}
    state::SP
end
const RLOp{SP<:ConcreteBasis} = Union{RaiseOp{SP}, LowerOp{SP}}

struct RaiseLowerOp{SP<:ConcreteBasis} <: AbstractOperator{Bases.Slater{SP}, Int}
    rlops::Vector{RLOp{SP}}

    RaiseLowerOp{SP}(rlops::Vector{RLOp{SP}}) where SP<:ConcreteBasis = new(rlops)
end
RaiseLowerOp{SP}(itr) where SP<:ConcreteBasis = RaiseLowerOp(collect(RLOp{SP}, itr))
RaiseLowerOp(rlops::Vector{RLOp{SP}}) where SP<:ConcreteBasis = RaiseLowerOp{SP}(rlops)

const AnyRLOp{SP<:ConcreteBasis} = Union{RaiseLowerOp{SP}, RaiseOp{SP}, LowerOp{SP}}

(op::RaiseOp)(args...) = applyop(op, args...)
(op::LowerOp)(args...) = applyop(op, args...)
(op::RaiseLowerOp)(args...) = applyop(op, args...)

ManyBody.reptype(O::Type{<:RLOp}) = innertype(basistype(O))
ManyBody.reptype(O::Type{<:RaiseLowerOp}) = Vector{RLOp{innertype(basistype(O))}}

function applyop!(op::RaiseOp{SP}, b::Bases.Slater{SP}) where SP<:ConcreteBasis
    sgn = create!(b, op.state)
    if sgn > 0
        b
    elseif sgn == 0
        States.Zero()
    else
        -b
    end
end
function applyop!(op::LowerOp{SP}, b::Bases.Slater{SP}) where SP<:ConcreteBasis
    sgn = annihil!(b, op.state)
    if sgn > 0
        b
    elseif sgn == 0
        States.Zero()
    else
        -b
    end
end
function applyop(op::RLOp{SP}, b::Bases.Slater{SP}) where SP<:ConcreteBasis
    ret = deepcopy(b)
    applyop!(op, ret)
end

applyop!(op::RaiseLowerOp{SP}, b::Bases.Slater{SP}) where SP<:ConcreteBasis =
    foldr(applyop!, op.rlops, init=b)
function applyop(op::RaiseLowerOp{SP}, b::Bases.Slater{SP}) where SP<:ConcreteBasis
    ret = deepcopy(b)
    applyop!(op, ret)
end

struct Raised{T}; val::T end
↑(x) = Raised(x)

A(p::ConcreteBasis) = LowerOp(p)
A(p::States.Bra{<:ConcreteBasis}) = RaiseOp(p.state)
A(p::Raised{<:ConcreteBasis}) = RaiseOp(p.val)

@generated function A(sps::Vararg{Union{T, States.Bra{T}, Raised{T}}}) where
                     {B<:ConcreteBasis, T<:Bases.Rep{B}}
    args = map(enumerate(sps)) do x
        i, U = x
        if U <: States.Bra
            :(RaiseOp(sps[$i].state))
        elseif U <: Raised
            if U.parameters[1] <: Bases.Wrapped{B}
                :(RaiseOp(inner(sps[$i].val)))
            else
                :(RaiseOp(sps[$i].val))
            end
        elseif U <: Bases.Wrapped{B}
            :(LowerOp(inner(sps[$i])))
        else
            :(LowerOp(sps[$i]))
        end
    end

    :(RaiseLowerOp{B}(RLOp{B}[$(args...)]))
end

refop(s::Bases.Wrapped) = refop(inner(s))
function refop(R::RefState{B}, s::Bases.Slater{B}) where B<:ConcreteBasis
    BI = indextype(B)
    RaiseLowerOp{B}(isocc(R, BI[I]) ? LowerOp{B}(B[I]) : RaiseOp{B}(B[I])
                    for I in eachindex(B)
                    if isocc(R, BI[I]) && ~s.bits[I] || ~isocc(R, BI[I]) && s.bits[I])
end

normord(a::RaiseLowerOp) = normord(RefStates.Vacuum(basistype(a)), a)
function normord(R::RefState{B}, a::RaiseLowerOp{B}) where B<:ConcreteBasis
    function comp(a, b)
        if a[2] && b[2]
            index(a[1]) >= index(b[1])
        elseif !a[2] && !b[2]
            index(a[1]) <= index(b[1])
        else
            a[2]
        end
    end

    xs = map(a.rlops) do b
        (b.state, if b isa RaiseOp
            ~isocc(R, b.state)
        else
            isocc(R, b.state)
        end)
    end

    p = sortperm(xs, lt=comp)
    sgn = levicivita(p)

    ret = RaiseLowerOp{B}(a.rlops[p])
    if sgn > 0
        ret
    elseif sgn < 0
#        ActionOperator{Bases.Slater{B}, Int}(x -> -ret(x))
        -ret
    else
        error("IMPOSSIBLE")
    end
end

Base.show(io::IO, x::RaiseOp) = print(io, "A($(x.state)')")
Base.show(io::IO, x::LowerOp) = print(io, "A($(x.state) )")
function Base.show(io::IO, x::RaiseLowerOp)
    f(a::RaiseOp) = "$(a.state)'"
    f(a::LowerOp) = "$(a.state) "

    print(io, "A(", f(x.rlops[1]))
    for op in x.rlops[2:end]
        print(io, ", $(f(op))")
    end
    print(io, ")")
end
function Base.show(io::IO, mime::MIME"text/plain", x::ArrayOperator)
    println(typeof(x))
    show(io, mime, x.rep)
end

Base.:(==)(a::ArrayOperator{B}, b::ArrayOperator{B}) where B<:ConcreteBasis = a.rep == b.rep

Base.copy(a::ArrayOperator) = typeof(a)(a.rep)

Base.zero(O::Type{<:ActionOperator}) = O(a -> States.Zero())
Base.zero(O::Type{<:FunctionOperator}) = O((a, b) -> zero(eltype(O)))
function Base.zero(O::Type{<:ArrayOperator})
    A = reptype(O)
    B = basistype(O)
    ds = innerdims(B)

    O(zero(similar(A, ds..., ds...)))
end
### No sensible zero for RLOp...
## And this is but vaguely sensible, and so shall be omitted
#Base.zero(O::RaiseLowerOp) = O*O

function Base.:+(a::ActionOperator{B}, b::ActionOperator{B}) where B<:ConcreteBasis
    Ta, Tb = eltype.((a, b))
    fa, fb = (a.rep, b.rep)
    ActionOperator{B, promote_type(Ta, Tb)}(x -> fa(x) + fb(x))
end

function Base.:+(a::FunctionOperator{B}, b::FunctionOperator{B}) where B<:ConcreteBasis
    Ta, Tb = eltype.((a, b))
    fa, fb = (a.rep, b.rep)
    FunctionOperator{B, promote_type(Ta, Tb)}((x, y) -> fa(x, y) + fb(x, y))
end

Base.:+(a::ArrayOperator{B}, b::ArrayOperator{B}) where B<:ConcreteBasis =
    ArrayOperator{B}(a.rep + b.rep)

Base.:+(a::AnyRLOp{B}, b::AnyRLOp{B}) where B<:ConcreteBasis =
    ActionOperator{Bases.Slater{B}, Int}(x -> a(x) + b(x))

Base.:-(a::ActionOperator) = (fa=a.rep; ActionOperator{basistype(a), eltype(a)}(x -> -fa(x)))
function Base.:-(a::FunctionOperator)
    fa = a.rep
    FunctionOperator{basistype(a), eltype(a)}((x, y) -> -fa(x, y))
end
Base.:-(a::ArrayOperator) = typeof(a)(-a.rep)
Base.:-(a::AnyRLOp) = ActionOperator{basistype(a), Int}(x -> -a(x))

Base.:-(a::AbstractOperator, b::AbstractOperator) = a + -b

@commutes function Base.:*(c::Number, a::ActionOperator)
    T = promote_type(typeof(c), eltype(a))
    f = a.rep
    ActionOperator{basistype(a), T}(x -> c*f(x))
end
@commutes function Base.:*(c::Number, a::FunctionOperator)
    T = promote_type(typeof(c), eltype(a))
    f = a.rep
    FunctionOperator{basistype(a), T,}((x, y) -> c*f(x, y))
end
@commutes Base.:*(c::Number, a::ArrayOperator) =
    ArrayOperator{basistype(a)}(c*a.rep)
@commutes Base.:*(c::Number, a::AnyRLOp) =
    ActionOperator{basistype(a), promote_type(typeof(c), Int)}(x -> c*a(x))


function Bases.widen(B::Type{<:ConcreteBasis}, op::ArrayOperator)
    SB = basistype(op)
    SBI = indextype(SB)
    ret = similar(reptype(op), innderdims(B)..., innerdims(B)...) |> zero

    for (I, J) in cartesian_pow(eachindex(SB), Val{2})
        I′, J′ = index(Bases.widen(B, SBI[I])), index(Bases.widen(B, SBI[J]))
        ret[I′, J′] = op.rep[I, J]
    end

    ArrayOperator{B}(ret)
end

function Bases.constrain(SB::Type{<:ConcreteBasis}, op::ArrayOperator)
    B = basistype(op)
    SBI = indextype(SB)
    ret = similar(reptype(op), innerdims(SB)..., innerdims(SB)...)

    for (I, J) in cartesian_pow(eachindex(SB), Val{2})
        I′, J′ = index(Bases.widen(B, SBI[I])), index(Bases.widen(B, SBI[J]))
        ret[I, J] = op.rep[I′, J′]
    end

    ArrayOperator{SB}(ret)
end

end # module Operators
