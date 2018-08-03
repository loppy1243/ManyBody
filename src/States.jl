@reexport module States
export State, VecState, CVecState, MaybeBasis

using ..Bases
using ..AbstractState
using Suppressor: @suppress_err

struct Zero <: AbstractState
struct State{B<:AbstractBasis, V<:AbstractVector} <: AbstractState
    vec::V
    State{B, V}(vec::V) where {B, V} = new(vec)
end
const VecState{B, T} = State{B, Vector{T}}
const CVecState{B} = VecState{B, Complex64}

function State{B, V}(x::B) where {B<:AbstractBasis, V<:AbstractVector}
    v = zero(V(dim(B)))
    v[index(x)] = oneunit(eltype(V))
    State{B, V}(v)
end
State{B}(x) where B = CVecState{B}(x)

dim(x::State) = length(x.vec)

const MaybeBasis{B<:AbstractBasis, V<:AbstractVector} = Union{B, State{B, V}}

# Add type arg
@suppress_err Base.vec(x::MaybeBasis{B, V}) where {B, V} = convert(State{B}, x).vec

Base.zero(::Type{State{B, V}}) where {B, V} = State{B, V}(zero(V(dim(B))))
Base.zero(x::State) = typeof(x)(zero(x.vec))

Base.convert(::Type{State{B, V1}}, x::State{B, V2}) where {B, V1, V2} =
    State{B, V1}(convert(V1, x.vec))
Base.convert(::Type{A}, x::B) where {B<:AbstractBasis, A<:State{B}} = A(x)
Base.convert(::Type{State{B}}, x::State{B}) where B = x
function Base.convert(::Type{B}, x::State{B}) where B
    T = eltype(x.vec)
    nzs = x.vec .!= zero(T)
    i = findfirst(nzs)
    cnt = count(zs)

    if cnt != 1 || x.vec[i] != oneunit(T)
        throw(InexactError())
    end

    B[i]
end

@suppress_err for op in (:+, :-), T in (:(MaybeBasis{B}), :(Bra{MaybeBasis{B}}))
    @eval begin
        function Base.$op(a::$T, b::$T) where B<:AbstractBasis
            a = convert(State{B}, a)
            b = convert(State{B}, b)
            ret = $op(a.vec, b.vec)

            State{B, typeof(ret)}(ret)
        end
    end
end

@suppress_err for T in (:(MaybeBasis{B}), :(Bra{MaybeBasis{B}}))
    @eval begin
        function Base.:*(x::Number, a::$T) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = a*vec(a)
        
            State{B, typeof(ret)}(ret)
        end
        function Base.:*(a::$T, x::Number) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = vec(a)*x
        
            State{B, typeof(ret)}(ret)
        end
        function Base.:/(a::$T, x::Number) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = vec(a)/x
        
            State{B, typeof(ret)}(ret)
        end
        function Base.:\(x::Number, a::$T) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = x \ vec(a)
        
            State{B, typeof(ret)}(ret)
        end
        Base.:+(a::$T) where B<:AbstractBasis = a

        function Base.:-(a::$T) where B<:AbstractBasis
            a = convert(State{B}, a)
            ret = -a.vec
            State{B, typeof(ret)}(ret)
        end
    end
end

struct Bra{S<:AbstractState}; state::S end

Base.:(==)(b1::Bra{S}, b2::Bra{S}) where S = b1.state == b2.state
Base.ctranspose(s::AbstractState) = Bra(s)
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
    a.vec'b.vec
end

end # module States
