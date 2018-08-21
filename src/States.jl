@reexport module States
export ArrayState, F64ArrayState, CF64ArrayState, VectorState, overlap
using ..Bases
using ..AbstractState, ..@commutes

abstract type Mixed{B<:ConcreteBasis, T} <: AbstractState end

struct Bra{S<:AbstractState} end
struct Zero <: AbstractState end

struct Scaled{B<:ConcreteBasis, T} <: Mixed{B, T}
    coeff::T
    state::B
end
const CF64Scaled{B<:ConcreteBasis} <: Scaled{B, ComplexF64}

Scaled(T::Type, b::ConcreteBasis) = Scaled{typeof(b), T}(one(T), b)
Scaled(c, b::ConcreteBasis) = Scaled{B, typeof(c)}(c, b)

struct ArrayState{B<:ConcreteBasis, T, A<:AbstractArray{T}} <: Mixed{B, T}
    coeffs::A

    function ArrayState{B, T, A}(coeffs::A) where
                       {B<:ConcreteBasis, T, A<:AbstractArray{T}}
        @assert size(coeffs) == innerdims(B)
        new(coeffs)
    end
end
const CF64ArrayState{B<:ConcreteBasis} =
    ArrayState{B, ComplexF64, Array{ComplexF64, N}}
const F64ArrayState{B<:ConcreteBasis} = ArrayState{B, Float64, Array{Float64, N}}

const CF64VectorState{B<:ConcreteBasis} = VectorState{B, ComplexF64, Vector{ComplexF64}}
const F64VectorState{B<:ConcreteBasis} = VectorState{B, Float64, Vector{FloatF64}}

ArrayState{B}(coeffs::AbstractArray) where B<:ConcreteBasis =
    ArrayState{B, eltype(coeffs), typeof(coeffs)}(coeffs)

function ArrayState(A::Type{<:AbstractArray}, b::ConcreteBasis)
    ret = zero(similar(A, innerdims(typeof(b))...))
    ret[index(b)] = oneunit(eltype(A))
    ArrayState{typeof(b), eltype(A), A}(ret)
end
ArrayState(A::Type{<:AbstractArray}, s::Scaled) = s.coeff*ArrayState(A, s.state)

Base.promote_rule(S1::Type{<:Scaled}, S2::Type{<:Scaled}) =
    Scaled{promote_type(basistype(S1), basistype(S2)), promote_type(eltype(S1), eltype(S2))}
function Base.promote_rule(S1::Type{<:ArrayState}, S2::Type{<:ArrayState})
    B = promote_type(basistype(S1), basistype(S2))
    A = typeof(similar(promote_type(reptype(S1), reptype(S2)), fill(0, rank(B))...))

    ArrayState{B, eltype(A), A}
end
Base.promote_rule(S::Type{<:Scaled{B}}, ::Type{B}) where B<:AbstractBasis = S
Base.promote_rule(S::Type{<:ArrayState{B}}, ::Type{B}) where B<:ConcreteBasis = S
##Base.promote_rule(S::Type{<:ArrayState}, B::Type{<:AbstractBasis}) =
function Base.promote_rule(S1::Type{<:ArrayState{B}}, S2::Type{<:Scaled{B}}) where
                 B<:ConcreteBasis
    T = promote_type(eltype(S1), eltype(S2))
    A = typeof(similar(reptype(S1), T, fill(0, rank(S1))))
    ArrayState{B, T, A}
end

Base.convert(S::Type{<:Scaled}, b) = S(b)
Base.convert(S::Type{<:ArrayState}, b) = S(b)
Base.convert(S1::Type{<:Scaled}, s2::Scaled) =
    S1(convert(eltype(S1), s2.coeff), convert(basistype(S1), s2.state))
function Base.convert(S1::Type{<:ArrayState}, s2::ArrayState)
    B1, B2 = basistype.((S1, s2))
    B1_I, B2_I = indextype.((B1, B2))
    A1 = reptype(S1)

    ret = similar(A1, innerdims(B1))
    for I in eachindex(B2)
        ret[index(B1_I[B2_I[I]])] = rep(s2)[I]
    end
end

#const Rep{B<:AbstractBasis, S<:Mixed{B}} = Union{Bases.Rep{B}, S}

Bases.dim(S::Type{<:Mixed}) = dim(basistype(S))
Bases.rank(S::Type{<:Mixed}) = rank(basistype(S))
ManyBody.basistype(::Type{<:Mixed{B}}) where B<:AbstractBasis = B

Base.eltype(::Type{<:Mixed{<:Any, T}}) where T = T

ManyBody.rep(s::Scaled) = (s.coeff, s.state)
ManyBody.rep(s::ArrayState) = s.coeffs
ManyBody.reptype(S::Type{<:Scaled}) = Tuple{eltype(S), basistype(S)}
ManyBody.reptype(S::Type{<:ArrayState{<:Any, <:Any, A}}) where A<:AbstractArray = A

Base.zero(s::Scaled) = typeof(s)(zero(eltype(s)), s.state)
Base.zero(s::ArrayState) = typeof(s)(zero(s.coeffs))
@generated Base.zero(S::Type{<:ArrayState}) =
    :(S(zero(similar(reptype(S), $(fill(dim(basistype(S)), rank(S))...)))))

### Update line
##############################################################################################

@commutes Base.:*(c::Number, b::ConcreteBasis) = Scaled(c, b)
@commutes Base.:*(c::Number, b::Bases.Neg) = Scaled(-c, inner(b))
@commutes Base.:*(c::Number, s::Scaled) = Scaled(c*s.coeff, s.state)
@commutes Base.:*(c::Number, s::ArrayState) = ArrayState{basistype(s)}(c*rep(s))

Base.:/(s::AbstractState, c::Number) = inv(c)*s
Base.:\(c::Number, s::AbstractState) = inv(c)*s

Base.:+(a::AbstractState) = a
Base.:+(a::AbstractState, b::AbstractState) = +(promote(a, b)...)

@commutes Base.:+(a::AbstractState, b::Bases.Neg) = a - inner(b)
Base.:+(a::B, b::B) where B<:ConcreteBasis = if a == b
    Scaled(2, a)
else
    ArrayState(Array{Int}, a) + ArrayState(Array{Int}, b)
end
    
Base.:+(a::S, b::S) where S<:Scaled = if a.state == b.state
    Scaled(a.coeff + b.coeff, a.state)
else
    A = Array{eltype(S)}
    ArrayState(A, a) + ArrayState(A, b)
end
Base.:+(a::ArrayState{B}, b::ArrayState{B}) where B<:ConcreteBasis =
    ArrayState{B}(rep(a) + rep(b))

Base.:-(a::Scaled) = Scaled(-a.coeff, a.state)
Base.:-(a::ArrayState) = ArrayState{basistype(a)}(-rep(a))
Base.:-(a::AbstractState, b::AbstractState) = a + -b

Base.adjoint(s::AbstractState) = Bra(s)
Base.:*(bra::Bra, ket::AbstractState) = overlap(bra.state, ket)

overlap(a::AbstractState, b::AbstractState) = overlap(promote(a, b)...)

overlap(a::B, b::B) where B<:AbstractBasis = a == b
overlap(a::Bases.MaybeNeg{B}, b::Bases.MaybeNeg{B}) where B<:ConcreteBasis =
    -(inner(a) == inner(b))
overlap(a::B, b::B) where B<:Bases.Neg = a == b
@commutes conj overlap(a::Bases.Neg{B}, b::MaybeMixed{B}) where B<:ConcreteBasis =
    -overlap(inner(a), b)

overlap(::Zero, ::Zero) = 0
@commutes overlap(::Zero, b::AbstractState) = 0
@commutes overlap(::Zero, b::Mixed) = zero(eltype(b))

overlap(a::S, b::S) where S<:Scaled = conj(a.coeff)*b.coeff*overlap(a.state, b.state)
overlap(a::ArrayState{B}, b::ArrayState{B}) where B<:ConcreteBasis =
    vec(rep(a))'vec(rep(b))

@commutes conj overlap(a::Scaled{B}, b::MaybeMixed{B}) where B<:ConcreteBasis =
    conj(a.coeff)*overlap(a.state, b)

@commutes conj overlap(a::B, b::ArrayState{B}) where B<:ConcreteBasis = rep(b)[index(a)]
@commutes conj overlap(a::Bases.Index{B}, b::ArrayState{B}) where B<:ConcreteBasis =
    rep(b)[index(a)]

end # module States
