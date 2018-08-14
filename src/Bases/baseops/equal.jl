Base.:(==)(x::B, y::B) where B<:AbstractBasis = index(x) == index(y)

Base.:(==)(a::AbstractArray, ::ZeroState) = all(iszero, a)
Base.:(==)(::ZeroState, b::AbstractArray) = all(iszero, b)

Base.:(==)(a::AbstractVector, b::AbstractBasis) = b == a
function Base.:(==)(a::AbstractBasis, b::AbstractVector)
    i = index(a)
    b[i] == oneunit(eltype(b)) && all(iszero(b[k]) for k in eachindex(b) if k != i)
end
