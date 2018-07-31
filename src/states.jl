struct State{B<:Basis, V<:AbstractVector{T}} <: State
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
    @eval Base.$op(a::S, b) where S<:State = S($op(a.vec, b))
    @eval Base.$op(a, b::S) where S<:State = S($op(a, b.vec))
    @eval Base.$op(a::S, b::S) where S<:State = S($op(a.vec, b.vec))
end
