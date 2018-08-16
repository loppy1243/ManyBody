abstract type MixedState{N, B<:AbstractBasis, T} <: AbstractState end

struct Zero <: AbstractState end

struct Scaled{B<:AbstractBasis, T} <: MixedState{1, B, T}
    coeff::T
    state::B
end
const CF64Scaled{B<:AbstractBasis} <: Scaled{B, ComplexF64}
Scaled{B, T}(b::B) where {B<:AbstractBasis, T} = Scaled{B, T}(one(T), b)
Scaled(::Type{T}, b::B) where {T, B<:AbstractBasis} = Scaled{B, T}(one(T), b)
Scaled(c, b::B) where B<:AbstractBasis = Scaled{B, typeof(c)}(c, b)
Scaled(::Type{T}, b::Neg) where T = Scaled(-one(T), inner(b))
Scaled(c, b::Neg) = Scaled(-c, inner(b))

struct ArrayState{N, B<:AbstractBasis, T, A<:AbstractArray{T, N}} <: MixedState{N, B, T}
    coeffs::A

    function ArrayState{N, B, T, A}(coeffs::A) where {N, B<:AbstractBasis, T, A<:AbstractArray{T, N}}
        @assert all(size(coeffs) .== dim(B))

        new(coeffs)
    end
end
const VectorState{B<:AbstractBasis, T, V<:AbstractVector{T}} = ArrayState{1, B, T, V}
const CF64VectorState{B<:AbstractBasis} = ArrayState{1, B, ComplexF64, Vector{ComplexF64}}
const CF64ArrayState{N, B<:AbstractBasis} = ArrayState{N, B, ComplexF64, Array{ComplexF64, N}}
const F64VectorState{B<:AbstractBasis} = ArrayState{1, B, Float64, Vector{Float64}}
const F64ArrayState{N, B<:AbstractBasis} = ArrayState{N, B, Float64, Array{Float64, N}}

ArrayState(B::Type{<:AbstractBasis}, coeffs::AbstractArray) =
    ArrayState{ndims(coeffs), B, eltype(coeffs), typeof(coeffs)}(coeffs)

function VectorState{B, V}(b::B) where {B<:AbstractBasis, V<:AbstractVector}
    ret = zero(similar(V, dim(B)))
    ret[index(b)] = oneunit(eltype(V))
    VectorState{B, V}(ret)
end

const Rep{B, S<:MixedState} = Union{Bases.Rep{B}, S}

dim(::Type{<:MixedState{N, B}}) where {N, B<:AbstractBasis} = dim(B)^N
rank(::Type{<:MixedState{N}}) where N = N

basistype(::Type{<:MixedState{<:Any, B}}) where B<:AbstractBasis = B
basisdim(S::Type{<:MixedState}) = dim(basistype(S))

Base.eltype(::Type{<:MixedState{<:Any, <:Any, T}}) where T = T

rep(s::Scaled) = (s.coeff, s.state)
rep(s::ArrayState) = s.coeffs
reptype(S::Type{<:Scaled}) = Tuple{eltype(S), basistype(S)}
reptype(S::Type{<:ArrayState{<:Any, <:Any, <:Any, A}}) where A<:AbstractArray = A

Base.zero(s::Scaled) = typeof(s)(zero(eltype(s)), s.state)
Base.zero(s::ArrayState) = typeof(s)(zero(s.coeffs))
@generated Base.zero(S::Type{<:ArrayState}) =
    :(S(zero(similar(reptype(S), $(fill(dim(basistype(S)), N)...)))))


Base.:*(c::Number, b::AbstractBasis) = Scaled(c, b)
Base.:*(b::AbstractBasis, c::Number) = Scaled(c, b)
Base.:/(b::AbstractBasis, c::Number) = Scaled(inv(c), b)
Base.:\(c::Number, b::AbstractBasis) = Scaled(inv(c), b)

Base.:+(a::Scaled{<:Rep{B}}, b::Rep{B}) where B<:AbstractBasis =
    a.state == b ? Scaled(a.coeff + one(eltype(a)), promote(a, b))
