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

for op in (:+, :-, :*, :/, :\)
    @eval function Base.$op(a::State{B}, b) where B
        ret = $op(a.vec, b)
        State{B, typeof(ret)}(ret)
    end
    @eval function Base.$op(a, b::State{B}) where B
        ret = $op(a, b.vec)
        State{B, typeof(ret)}(ret)
    end
    @eval function Base.$op(a::S1, b::S2) where {B, S1<:State{B}, S2<:State{B}}
        ret = $op(a.vec, b.vec)
        State{B, typeof(ret)}(ret)
    end
end

for op in (:+, :-)
    @eval Base.$op(a::State) = typeof(a)($op(a.vec))
end

overlap(a::S, b::S) where {B, S<:State{B}} = a.vec'b.vec
