Base.:(==)(x::AbstractBasis, y::AbstractBasis) = false
Base.:(==)(x::B, y::B) where B<:AbstractBasis = index(x) == index(y)
Base.:(==)(a::AbstractArray, ::ZeroState) = all(iszero, a)
Base.:(==)(::ZeroState, b::AbstractArray) = all(iszero, b)
function Base.:(==)(a::AbstractBasis, b::AbstractVector)
    i = index(a)
    b[i] == oneunit(eltype(b)) && all(iszero(b[k]) for k in eachindex(b) if k != i)
end
Base.:(==)(a::AbstractVector, b::AbstractBasis) = b == a

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

Base.adjoint(b::AbstractBasis) = Bra(b)
Base.adjoint(b::Bra) = b.state

Base.:*(a::Bra{B}, b::B) where B<:AbstractBasis = a == b
Base.:*(::Bra{ZeroState}, ::Union{<:AbstractVector, <:AbstractBasis}) = 0
Base.:*(::Union{<:Adjoint{<:Any, <:AbstractVector}, Bra{<:AbstractBasis}}, ::ZeroState) = 0
Base.:*(::Bra{ZeroState}, ::ZeroState) = 0
Base.:*(a::Bra{<:AbstractBasis}, b::AbstractVector) = b[index(a)]
Base.:*(a::Adjoint{<:Any, <:AbstractVector}, b::AbstractBasis) = conj(b[index(a)])

for op in (:*, :/)
    @eval Base.$op(a::Number, b::AbstractBasis) = $op(a, Vector{typeof(a)}(b))
    @eval Base.$op(::Number, ::ZeroState) = ZeroState()
end
for op in (:*, :\)
    @eval Base.$op(a::AbstractBasis, b::Number) = $op(Vector{typeof(a)}(a), b)
    @eval Base.$op(::ZeroState, ::Number) = ZeroState()
end

for op in (:+, :-); @eval begin
    Base.$op(a::Union{<:AbstractArray, <:AbstractBasis}, ::ZeroState) = a
    Base.$op(::ZeroState, b::Union{<:AbstractArray, <:AbstractBasis}) = $op(b)
    Base.$op(::ZeroState, ::ZeroState) = ZeroState()

    Base.$op(a::B, b::B) where B<:AbstractBasis = $op(Vector(a), Vector(b))
    Base.$op(a::AbstractVector, b::AbstractBasis) = $op(a, convert(typeof(a), b))
    Base.$op(a::AbstractBasis, b::AbstractVector) = $op(convert(typeof(b), a), b)
    Base.$op(a::B, b::Sub{B}) where B<:AbstractBasis = $op(Vector(a), Vector(convert(B, b)))
    Base.$op(a::Sub{B}, b::B) where B<:AbstractBasis = $op(Vector(convert(B, a)), Vector(b))
end; end
Base.:+(a::AbstractBasis) = a
Base.:-(a::AbstractBasis) = -Vector(a)

Base.promote_rule(::Type{V}, ::Type{<:AbstractBasis}) where {V<:AbstractVector} = V

## OMG this worked. Please add actual tests.
@generated function Base.:*(xs::Vararg{Union{AbstractArray, AbstractBasis}})
    numdims = sum(xs) do x
        if x <: AbstractArray
            ndims(x)
        else
            1
        end
    end
    dim_exprs = map(enumerate(xs)) do X
        k, x = X
        if x <: AbstractArray
            Expr(:..., :(size(xs[$k])))
        else
            :(dim(typeof(xs[$k])))
        end
    end

    i = 1
    j = 1
    exprs = map(enumerate(xs)) do X
        k, x = X
        if x <: AbstractArray
            ret = :(xs[k][CartesianIndex(J[$(i:i+ndims(x)-1)])])
            i += ndims(x)
            ret
        else
            ret = :(J[$i] == ixs[$j])
            j += 1
            i += 1
            ret
        end
    end

    ixs_expr = [:(index(xs[$k])) for k in 1:lastindex(xs) if xs[k] <: AbstractBasis]

    quote
        ixs = [$(ixs_expr...)]
        ret = zeros($(dim_exprs...))
        for I in CartesianIndices(ret)
            J = Tuple(I)
            ret[I] = *($(exprs...))
        end
        ret
    end
end
