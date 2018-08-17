Base.adjoint(b::AbstractBasis) = Bra(b)
Base.adjoint(b::Bra) = b.state

#Base.:*(a::Bra{B}, b::B) where B<:AbstractBasis = a == b
Base.:*(a::Bra{<:Rep{B}}, b::Rep{B}) where B<:AbstractBasis = inner(a.state) == inner(b)
#Base.:*(::Bra{ZeroState}, ::Union{<:AbstractVector, <:AbstractBasis}) = 0
#Base.:*(::Union{<:Adjoint{<:Any, <:AbstractVector}, Bra{<:AbstractBasis}}, ::ZeroState) = 0
Base.:*(::Bra{ZeroState}, ::ZeroState) = 0
#Base.:*(a::Bra{<:AbstractBasis}, b::AbstractVector) = b[index(a.state)]
#Base.:*(a::Adjoint{<:Any, <:AbstractVector}, b::AbstractBasis) = conj(a[index(b)])

for op in (:*, :/); @eval begin
#    Base.$op(a::Number, b::AbstractBasis) = $op(a, Vector{typeof(a)}(b))
    Base.$op(::Number, ::ZeroState) = ZeroState()
end; end
for op in (:*, :\); @eval begin
#   Base.$op(a::AbstractBasis, b::Number) = $op(Vector{typeof(a)}(a), b)
    Base.$op(::ZeroState, ::Number) = ZeroState()
end; end

for op in (:+, :-); @eval begin
#   Base.$op(a::Union{<:AbstractArray, <:AbstractBasis}, ::ZeroState) = a
#   Base.$op(::ZeroState, b::Union{<:AbstractArray, <:AbstractBasis}) = $op(b)
    Base.$op(::ZeroState, ::ZeroState) = ZeroState()

#    Base.$op(a::B, b::B) where B<:AbstractBasis = $op(Vector(a), Vector(b))
#    Base.$op(a::AbstractVector, b::AbstractBasis) = $op(a, convert(typeof(a), b))
#    Base.$op(a::AbstractBasis, b::AbstractVector) = $op(convert(typeof(b), a), b)
#    Base.$op(a::B, b::Sub{B}) where B<:AbstractBasis = $op(Vector(a), Vector(convert(B, b)))
#    Base.$op(a::Sub{B}, b::B) where B<:AbstractBasis = $op(Vector(convert(B, a)), Vector(b))
end; end
Base.:+(a::AbstractState) = a
Base.:-(::ZeroState) = ZeroState()
#Base.:-(a::AbstractBasis) = -Vector(a)

#if !hasmethod(Base.:*, (Vararg{AbstractArray},))
#   Base.:*(xs::Vararg{AbstractArray}) = MethodError(Base.:*, (xs...,)) |> throw
#end
## OMG this worked. Please add actual tests.
#generated function Base.:*(xs::Vararg{Union{AbstractArray, AbstractBasis}})
#    numdims = sum(xs) do x
#       if x <: AbstractArray
#            ndims(x)
#        else
#            1
#        end
#    end
#    dim_exprs = map(enumerate(xs)) do X
#        k, x = X
#       if x <: AbstractArray
#            Expr(:..., :(size(xs[$k])))
#        else
#            :(dim(typeof(xs[$k])))
#        end
#    end
#
#    i = 1
#    j = 1
#    exprs = map(enumerate(xs)) do X
#        k, x = X
#       if x <: AbstractArray
#            ret = :(xs[k][CartesianIndex(J[$(i:i+ndims(x)-1)])])
#            i += ndims(x)
#            ret
#        else
#            ret = :(J[$i] == ixs[$j])
#            j += 1
#            i += 1
#            ret
#        end
#    end
#
#    ixs_expr = [:(index(xs[$k])) for k in 1:lastindex(xs) if xs[k] <: AbstractBasis]
#
#    quote
#        ixs = [$(ixs_expr...)]
#        ret = zeros($(dim_exprs...))
#        for I in CartesianIndices(ret)
#            J = Tuple(I)
#            ret[I] = *($(exprs...))
#        end
#        ret
#    end
#end
