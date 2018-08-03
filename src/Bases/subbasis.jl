export subindices, subindexes
using Combinatorics: combinations
using Loppy.Util: cartesian_pow
using ..SpinMod

struct Sub{B<:AbstractBasis, T} <: AbstractBasis
    state::B
end
const MaybeSub{B<:AbstractBasis} = Union{B, <:Sub{B}}

Base.convert(::Type{B}, x::Sub{B}) where B = x.state
Base.convert(::Type{C}, x::Sub{B}) where {C, B<:C} = x.state
Base.convert(::Type{SB}, x::B) where {B, SB<:Sub{B}} = SB(x)

Base.:(==)(x::SB, y::SB) where SB<:Sub = x.state == y.state

dim(::Type{SB}) where SB<:Sub = length(subindices(SB))
index(x::Sub) = findfirst(basis(typeof(x)), x)
indexbasis(::Type{SB}, i) where {B, SB<:Sub{B}} = basis(SB)[i]

subindices(::Type{B}) where B<:AbstractBasis = indices(B)
subindexes(::Type{SB}) where SB<:Sub = map(Index{SB}, subindices(SB))

macro defSubBasis(ty_expr::Expr, expr)
    @assert ty_expr.head == :(<:)

    name(x::Symbol) = x
    name(x::Expr) = name(x.args[1])

    new_ty, base_ty = ty_expr.args[1:2]
    tys = if new_ty.head == :curly
        new_ty.args[2:end]
    else
        nothing
    end
    tys_esc = map(esc, tys)
    add_where(x) = tys == nothing ? x : Expr(:where, x, tys_esc...)
        
    subbasis_ty_esc = esc(gensym(string(name(new_ty))))

    new_ty_esc, base_ty_esc = esc.((new_ty, base_ty))

    subindices_expr(x) = add_where(:(Bases.subindices($x::Type{$new_ty_esc})))
    basis_expr = add_where(:(Bases.basis(::Type{$new_ty_esc})))

    quote
        struct $subbasis_ty_esc end
        const $new_ty_esc = Sub{$base_ty_esc, $subbasis_ty_esc}

        @generated $basis_expr = begin
            x = (() -> $(esc(expr)))()
            x = convert(Vector{$base_ty_esc}, x)
            @assert allunique(x)
            map($new_ty_esc, x)
        end
        @generated $(subindices_expr(:T)) = map(index, basis(T.parameters[1]))
    end
end
