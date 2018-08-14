Base.convert(::Type{B}, s::Sub{Index{B}}) where B<:AbstractBasis = convert(B, s.state)
Base.convert(::Type{B}, s::Index{SB}) where {B<:AbstractBasis, SB<:Sub{B}} = convert(B, SB[index(s)])
function Base.convert(::Type{V}, s::AbstractBasis) where V<:AbstractVector
    ret = V(undef, dim(typeof(s)))
    ret .= zero(eltype(V))
    ret[index(s)] = oneunit(eltype(V))

    ret
end
function Base.convert(::Type{B}, v::AbstractVector) where B<:AbstractBasis
    nzs = v .!= zero(eltype(v))
    i = findfirst(!iszero, nzs)
    if count(nzs) != 1 || v[i] != oneunit(eltype(v))
        InexactError() |> throw
    end

    B[i]
end

Base.Vector{T}(b::AbstractBasis) where T = convert(Vector{T}, b)
Base.Vector(b::AbstractBasis) where T = convert(Vector{ComplexF64}, b)

Base.promote_rule(::Type{V}, ::Type{<:AbstractBasis}) where {V<:AbstractVector} = V
