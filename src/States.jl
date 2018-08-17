module States
using ..Bases
using ..AbstractState

abstract type MixedState{N, B<:AbstractBasis, T} <: AbstractState end

struct Bra{S<:AbstractState} end
struct Zero <: AbstractState end

struct Scaled{B<:AbstractBasis, T} <: MixedState{B, T}
    coeff::T
    state::B

    function Scaled{B, T}(c::T, b::B) where {B<:AbstractBasis, T}
        @assert !(B <: Bases.Neg || B <: Bases.Product{1})
        new(c, b)
    end
end
const CF64Scaled{B<:AbstractBasis} <: Scaled{B, ComplexF64}

Scaled(T::Type, b::AbstractBasis) = Scaled{typeof(b), T}(one(T), b)
Scaled(c, b::AbstractBasis) = Scaled{B, typeof(c)}(c, b)
Scaled(T::Type, b::Bases.Neg{<:AbstractBasis}) = Scaled(-one(T), inner(b))
Scaled(c, b::Bases.Neg{<:AbstractBasis}) = Scaled(-c, inner(b))
Scaled(T::Type, b::Bases.Product{1}) = Scaled(T, inner(b))
Scaled(c, b::Bases.Product{1}) = Scaled(c, inner(b))

## Allow B<:Bases.Rep{<:Bases.Product{N}}?
struct ArrayState{N, B<:Bases.Product{N}, T, A<:AbstractArray{T, N}} <: MixedState{B, T}
    coeffs::A

    function ArrayState{N, P, T, A}(coeffs::A) where
                                   {N, P<:Bases.Product{N}, T, A<:AbstractArray{T, N}}
        @assert size(coeffs) == innerdims(P)

        new(coeffs)
    end
end
const VectorState{B<:AbstractBasis, T, V<:AbstractVector{T}} =
    ArrayState{1, Bases.Product{1, B}, T, V}

const CF64ArrayState{N, B<:Bases.Product{N}} =
    ArrayState{N, B, ComplexF64, Array{ComplexF64, N}}
const F64ArrayState{N, B<:Bases.Product} = ArrayState{N, B, Float64, Array{Float64, N}}

const CF64VectorState{B<:AbstractBasis} = VectorState{B, ComplexF64, Vector{ComplexF64}}
const F64VectorState{B<:AbstractBasis} = VectorState{B, Float64, Vector{FloatF64}}

ArrayState(B::Type{<:Bases.Product}, coeffs::AbstractArray) =
    ArrayState{ndims(coeffs), B, eltype(coeffs), typeof(coeffs)}(coeffs)
ArrayState(B::Type{<:AbstractBasis}, coeffs::AbstractVector) =
    ArrayState(Bases.Product{1, Tuple{B}}, coeffs)

ArrayState(A::Type{<:AbstractArray}, bs::Vararg{AbstractBasis}) =
    ArrayState(A, Bases.Product{length(bs), typeof(bs)}(bs))
function ArrayState(A::Type{<:AbstractArray}, b::Bases.Product)
    ret = zero(similar(A, innerdims(B)...))
    ret[CartesianIndex(innerindices(b))] = oneunit(eltype(A))
    ArrayState{rank(b), typeof(b), eltype(A), A}(ret)
end
ArrayState(V::Type{<:AbstractVector}, b::AbstractBasis) = ArrayState(V, Bases.Product(b))
ArrayState(A::Type{<:AbstractArray}, s::Scaled) = s.coeff*ArrayState(A, s.state)

Base.promote_rule(S1::Type{<:Scaled}, S2::Type{<:Scaled}) =
    Scaled{promote_type(basistype(S1), basistype(S2)), promote_type(eltype(S1), eltype(S2))}
function Base.promote_rule(S1::Type{<:ArrayState{N, P}}, S2::Type{<:ArrayState{N, P}}) where
                          {N, P<:Bases.Product{N}}
    A = promote_type(reptype(S1), reptype(S2))
    ArrayState{N, P, eltype(A), A}
end
#function Base.promote_rule(S1::Type{<:ArrayState{N}}, S2::Type{<:ArrayState{N}}) where N
#    B = promote_type(basistype(S1), basistype(S2))
#    A = promote_type(reptype(S1), reptype(S2))
#    ArrayState{N, B, eltype(A), A}
#end
Base.promote_rule(S::Type{<:Scaled{B}}, ::Type{B}) where B<:AbstractBasis = S
Base.promote_rule(S::Type{<:ArrayState{1, B}}, ::Type{<:Scaled{B}}) where B<:AbstractBasis = S
Base.promote_rule(S::Type{<:ArrayState{N, P}}, ::Type{<:Scaled{P}}) where
                 {N, P<:Bases.Product{N}} = S

Base.convert(S::Type{<:Scaled}, b) = S(b)
Base.convert(S::Type{<:ArrayState}, b) = S(b)
Base.convert(S1::Type{Scaled{<:Bases.Rep{B}}}, s2::Scaled{<:Bases.Rep{B}}) where
            B<:AbstractBasis =
    S1(convert(eltype(S1), s2.coeff), s2.state)
Base.convert(S1::Type{<:ArrayState{N, P}}, s2::ArrayState{N, P}) where {N, P<:Bases.Product{N}} =
    S1(convert(reptype(S1), rep(s2)))
#Base.convert(S1::Type{<:Scaled}, s2::ArrayState) = ...

const Rep{B<:AbstractBasis, S<:MixedState{B}} = Union{Bases.Rep{B}, S}

dim(S::Type{<:MixedState}) = dim(basistype(S))
rank(S::Type{<:MixedState}) = rank(S)

basistype(::Type{<:MixedState{<:Any, B}}) where B<:AbstractBasis = B
basisdim(S::Type{<:MixedState}) = dim(basistype(S))

Base.eltype(::Type{<:MixedState{<:Any, T}}) where T = T

rep(s::Scaled) = (s.coeff, s.state)
rep(s::ArrayState) = s.coeffs
reptype(S::Type{<:Scaled}) = Tuple{eltype(S), basistype(S)}
reptype(S::Type{<:ArrayState{<:Any, <:Any, <:Any, A}}) where A<:AbstractArray = A

Base.zero(s::Scaled) = typeof(s)(zero(eltype(s)), s.state)
Base.zero(s::ArrayState) = typeof(s)(zero(s.coeffs))
@generated Base.zero(S::Type{<:ArrayState}) =
    :(S(zero(similar(reptype(S), $(fill(dim(basistype(S)), N)...)))))

macro commutes(expr::Expr)
    _commutes(:identity, expr)
end
macro commute(f, expr::Expr)
    _commutes(f, expr)
end
function _commutes(f, expr::Expr)
    # Just assume we have a function expr

    commed_expr = deepcopy(expr)

    if commed_expr.args[1].head === :where
        reverse!(@view commed_expr.args[1].args[1].args[2:end])
    else
        reverse!(@view commed_expr.args[1].args[2:end])
    end

    if f !== :identity
        commed_expr.args[2] = :($f(commed_expr.args[2]))
    end

    quote
        $(esc(expr))
        $(esc(commed_expr))
    end
end

@commutes Base.:*(c::Number, b::AbstractBasis) = Scaled(c, b)
@commutes Base.:*(c::Number, s::Scaled) = Scaled(c*s.coeff, s.state)
@commutes Base.:*(c::Number, s::ArrayState) = ArrayState(basistype(s), c*rep(s))

Base.:/(s::AbstractState, c::Number) = inv(c)*s
Base.:\(c::Number, s::AbstractState) = inv(c)*s

Base.:+(a::AbstractState) = a
Base.:+(a::AbstractState, b::AbstractState) = +(promote(a, b)...)
Base.:+(a::B, b::B) where B<:AbstractBasis = if a == b
    Scaled(2, a)
else
    ArrayState(Array{Int}, a) + ArrayState(Array{Int}, b)
end
    
Base.:+(a::Scaled{B, T}, b::Scaled{B, T}) where {B<:AbstractBasis, T} = if a.state == b.state
    Scaled(a.coeff + b.coeff, a.state)
else
    ArrayState(Array{T}, a) + ArrayState(Array{T}, b)
end
Base.:+(a::ArrayState{N, P}, b::ArrayState{N, P}) where {N, P<:Bases.Product{N}} =
    ArrayState(P, rep(a) + rep(b))

Base.:-(a::Scaled) = Scaled(-a.coeff, a.state)
Base.:-(a::ArrayState) = ArrayState(basistype(a), -rep(a))
Base.:-(a::AbstractState, b::AbstractState) = a + -b

Base.adjoint(s::AbstractState) = Bra(s)
Base.:*(bra::Bra, ket::AbstractState) = overlap(bra.state, ket)

overlap(a::MaybeMixed{B}, b::MaybeMixed{B}) where B<:AbstractBasis = overlap(promote(a, b)...)

overlap(::Zero, ::Zero) = 0
overlap(a::B, b::B) where B<:AbstractBasis = a == b
overlap(a::Bases.Rep{B}, b::Bases.Rep{B}) where B<:AbstractBasis = inner(a) == inner(b)
overlap(a::Bases.MaybeNeg{B}, b::Bases.MaybeNeg{B}) where B<:AbstractBasis = -(inner(a) == inner(b))
@commutes conj overlap(a::Bases.Neg{B}, b::MaybeMixed{B}) where B<:AbstractBasis = -overlap(inner(a), b)

@commutes overlap(::Zero, b::AbstractState) = 0
@commutes overlap(::Zero, b::MixedState) = zero(eltype(b))

@commutes conj overlap(a::Scaled{<:Bases.Rep{B}}, b::Scaled{<:Bases.Rep{B}}) =
    conj(a.coeff)*b.coeff*overlap(a.state, b.state)
@commutes conj overlap(a::ArrayState{N, P}, b::ArrayState{N, P}) where {N, P<:Bases.Product{N}} =
    vec(rep(a))'vec(rep(b))

@commutes conj overlap(a::Scaled{<:Bases.Rep{B}}, b::MaybeMixed{B}) where B<:AbstractBasis =
    conj(a.coeff)*overlap(a.state, b)

@commutes conj overlap(a::Bases.Rep{B}, b::VectorState{B}) where B<:AbstractBasis =
    rep(b)[index(inner(a))]
@commutes conj overlap(a::Bases.Index{B}, b::VectorState{B}) where B<:AbstractBasis =
    rep(b)[index(a)]
## Make this more efficient in the Bases.Index{P} case somehow?
@commutes conj overlap(a::Bases.Rep{P}, b::ArrayState{N, P}) where {N, P<:Bases.Product{N}} =
    rep(b)[CartesianIndex(innerindices(inner(a.state)))]

end # module States
