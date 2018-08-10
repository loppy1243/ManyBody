@reexport module States
export State, VecState, CVecState, MaybeBasis, Bra

using ..Bases
using ..AbstractState
using Suppressor: @suppress_err

struct Bra{S<:AbstractState}; state::S end
struct Zero <: AbstractState end
struct State{B<:AbstractBasis, A<:AbstractArray} <: AbstractState
    coeffs::A
    State{B, A}(coeffs::A) where {B, A} = new(coeffs)
end
const ArrayState{B<:AbstractBasis, T, N} = State{B, Array{T, N}}
const CArrayState{B<:AbstractBasis, N} = State{B, Array{ComplexF64, N}}
const VecState{B<:AbstractBasis, T} = State{B, Array{T, 1}}
const CVecState{B<:AbstractBasis} = State{B, Array{ComplexF64, 1}}

const ZERO = Zero()

function State{B, A}(x::B) where {B<:AbstractBasis, A<:AbstractArray}
    coeffs = zero(A(undef, fill(dim(B), ndims(A))))
    coeffs[CartesianIndex(index(x))] = oneunit(eltype(A))
    State{B, A}(coeffs)
end

dim(x::State) = length(x.coeffs)
rank(x::State) = ndims(x.coeffs)

const MaybeBasis{B<:AbstractBasis, A<:AbstractArray} = Union{B, State{B, A}}

# Add type arg
@suppress_err Base.vec(x::MaybeBasis{B, A}) where {B, A<:AbstractVector} = convert(State{B}, x).vec

Base.zero(::Type{State{B, V}}) where {B, V} = State{B, V}(zero(V(dim(B))))
Base.zero(x::State) = typeof(x)(zero(x.vec))

Base.convert(::Type{State{B, A1}}, x::State{B, A2}) where
            {B, T, A1<:AbstractArray{<:Any, T}, A2<:AbstractArray{<:Any, T}} =
    State{B, A1}(convert(A1, x.coeffs))
Base.convert(::Type{A}, x::B) where {B<:AbstractBasis, A<:State{B}} = A(x)
#Base.convert(::Type{State{B}}, x::State{B}) where B = x
function Base.convert(::Type{B}, x::State{B, V}) where {B<:AbstractBasis, V<:AbstractVector}
    T = eltype(x.coeffs)
    nzs = x.coeffs .!= zero(T)
    i = findfirst(!iszero, nzs)
    cnt = count(nzs)

    if cnt != 1 || x.coeffs[i] != oneunit(T)
        throw(InexactError())
    end

    B[i]
end
function Base.convert(::Type{NTuple{N, B}}, x::State{B, A}) where
                     {B<:AbstractBasis, A<:AbstractArray{<:Any, N}}
    T = eltype(x.coeffs)
    nzs = x.coeffs .!= zero(T)
    ix = findfirst(!iszero, nzs)
    cnt = count(nzs)

    if cnt != 1 || x.coeffs[ix] != oneunit(T)
        throw(InexactError())
    end

    (map(i -> B[i], convert(Tuple, ix))...)
end

Base.promote_rule(::Type{SS}, ::Type{SB}) where
                 {B<:AbstractBasis, SS<:State{<:Bases.Sub{B}}, SB<:State{B}} = SB
Base.promote_rule(::Type{<:Bases.MaybeSub{B}}, ::Type{S}) where
                 {B<:AbstractBasis, S<:State{B}} = S

@suppress_err for op in (:+, :-), T in (:(MaybeBasis{B}), :(Bra{MaybeBasis{B}}))
    @eval begin
        function Base.$op(a::$T, b::$T) where B<:AbstractBasis
            a = convert(State{B}, a)
            b = convert(State{B}, b)
            ret = $op(a.coeffs, b.coeffs)

            State{B, typeof(ret)}(ret)
        end
    end
end

@suppress_err for T in (:(MaybeBasis{B}), :(Bra{<:MaybeBasis{B}}))
    @eval begin
        function Base.:*(x::Number, a::$T) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = a*a.coeffs
        
            State{B, typeof(ret)}(ret)
        end
        function Base.:*(a::$T, x::Number) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = coeffs(a)*x
        
            State{B, typeof(ret)}(ret)
        end
        function Base.:/(a::$T, x::Number) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = a.coeffs/x
        
            State{B, typeof(ret)}(ret)
        end
        function Base.:\(x::Number, a::$T) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = x \ a.coeffs
        
            State{B, typeof(ret)}(ret)
        end
        Base.:+(a::$T) where B<:AbstractBasis = a

        function Base.:-(a::$T) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = -a.coeffs
            State{B, typeof(ret)}(ret)
        end
    end
end

Base.:(==)(b1::Bra{S}, b2::Bra{S}) where S = b1.state == b2.state
Base.adjoint(s::AbstractState) = Bra(s)
Base.:*(bra::Bra, ket::AbstractState) = overlap(bra.state, ket)
Base.:*(x::Number, a::Zero) = Zero()
Base.:*(a::Zero, x::Number) = Zero()

overlap(a::B, b::B) where B<:AbstractBasis = a == b
overlap(a::Zero, b::S) where S<:AbstractState = 0
overlap(a::S, b::Zero) where S<:AbstractState = 0
overlap(a::Zero, b::Zero) = 0
overlap(a::MaybeBasis{B}, b::Bases.Sub{B}) where B<:AbstractBasis = overlap(a, convert(B, b))
overlap(a::Bases.Sub{B}, b::MaybeBasis{B}) where B<:AbstractBasis = overlap(convert(B, a), b)
function overlap(a::MaybeBasis{B}, b::MaybeBasis{B}) where B<:AbstractBasis
    a = convert(State{B}, a)
    a.coeffs'b.coeffs
end

end # module States
