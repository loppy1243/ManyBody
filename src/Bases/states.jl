export State, VecState, CVecState

struct State{B<:Basis, V<:AbstractVector{T}} <: AbstractState
    vec::V
end
const VecState{B, T} = State{B, Vector{T}}
const CVecState{B} = VecState{B, Complex64}

function State{B, V}(x::B) where {B, V}
    v = zeros(dim(B))
    v[index(x)] = oneunit(eltype(V))
    State{B, V}(v)
end
State{V}(x::B) where {B, V} = State{B, V}(x)

dim(x::State) = length(x.vec)

Base.vec(x::BasisOrState{B, V}) where {B, V} = convert(State{B, V}, x).vec

Base.zero(::Type{State{B, V}}) where {B, V} = State{B, V}(zero(V(dim(B))))
Base.zero(x::State) = typeof(x)(zero(x.vec))

Base.convert(::Type{State{B, V1}}, x::State{B, V2}) = State{B, V1}(convert(V1, x.vec))
Base.convert(::Type{A}, x::B) where {B, A<:State{B}} = A(x)
function Base.convert(::Type{B}, x::State{B})
    T = eltype(x.vec)
    nzs = x.vec .!= zero(T)
    i = findfirst(nzs)
    cnt = count(zs)

    if cnt != 1 || x.vec[i] != oneunit(T)
        throw(InexactError())
    end

    B[i]
end

const BasisOrState{B, V} = Union{B, State{B, V}}
for op in (:+, :-, :*, :/, :\)
    @eval begin
        function Base.$op(a::BasisOrState{B, V}, b) where {B, V}
            a = convert(State{B, V}, a)
            ret = $op(a.vec, b)
            State{B, typeof(ret)}(ret)
        end
        function Base.$op(a, b::BasisOrState{B, V}) where {B, V}
            b = convert(State{B, V}, b)
            ret = $op(a, b.vec)
            State{B, typeof(ret)}(ret)
        end
        function Base.$op(a::BasisOrState{B, V1}, b::BasisOrState{B, V2}) where {B, V2, V2}
            a = convert(State{B, V1}, a)
            b = convert(State{B, V2}, b)
            ret = $op(a.vec, b.vec)
            State{B, typeof(ret)}(ret)
        end
    end
end

Base.:+(a::BasisOrState) = a
for op in (:-, :conj)
    @eval function Base.$op(a::BasisOrState{B, V}) where {B, V}
        a = convert(State{B, V}, a)
        ret = $op(a.vec)
        State{B, typeof(ret)}(ret)
    end
end

overlap(a::S, b::S) where {B, S<:State{B}} = a.vec'b.vec
