using ..SpinMod

struct Sub{B<:AbstractBasis, T} <: AbstractBasis
    state::B
end
const MaybeSub{B<:AbstractBasis} = Union{B, Sub{B}}

#Base.convert(::Type{B}, x::Sub{B}) where B<:AbstractBasis = x.state
Base.convert(::Type{C}, x::Sub{B}) where {C<:AbstractBasis, B<:C} = x.state
Base.convert(::Type{SB}, x::B) where {B, SB<:Sub{B}} = SB(x)

Base.:(==)(x::S, y::S) where S<:Sub = inner(x) == inner(y)
Base.:(==)(x::Sub{B}, y::Sub{B}) where B<:AbstractBasis = false
Base.:(==)(x::MaybeSub{B}, y::MaybeSub{B}) where B<:AbstractBasis = inner(x) == inner(y)
Base.promote_rule(::Type{<:Sub{B}}, ::Type{B}) where B<:AbstractBasis = B

innertype(::Type{<:Sub{B}}) where B<:AbstractBasis = B
inner(s::Sub{<:AbstractBasis}) = s.state

get_tys(s::Symbol) = [s]
get_tys(x::Expr) = if x.head == :curly
    reduce(vcat, map(get_tys, x.args[2:end]))
else
    []
end
function get_free_tys(tys, x::Expr)
    @assert x.head === :curly
    ret = []
    for y in x.args[2:end]
        @assert y isa Symbol || y isa Expr && y.head === :(<:)
        if !(y in tys)
            push!(ret, y)
        end
    end

    ret
end
macro defSub(ty_expr::Expr, expr)
    _defSub(:Provided, ty_expr, expr)
end
## TODO: Currently assumes basis is Provided()
#macro defSub(gen, ty_expr::Expr, expr)
#    _defSub(gen, ty_expr, expr)
#end
function _defSub(gen, ty_expr, expr)
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
        
    free_tys = get_free_tys(get_tys(base_ty), new_ty)
    subbasis_ty_esc = esc(Expr(:curly, gensym(string(name(new_ty))), free_tys...))

    new_ty_esc, base_ty_esc = esc.((new_ty, base_ty))

    basis_expr = add_where(:(Bases.basis(::Type{$new_ty_esc})))
    generation_expr = add_where(:(Bases.Generation(::Type{$new_ty_esc})))
    dict_key_ty = add_where(:(Type{$new_ty_esc}))
    dict_val_ty = add_where(:(Vector{$new_ty_esc}))

    quote
        struct $subbasis_ty_esc end
        const $new_ty_esc = Sub{$base_ty_esc, $subbasis_ty_esc}

        $generation_expr = Bases.$gen()

        let _BASIS_CACHE = Dict{$dict_key_ty, $dict_val_ty}()
            $basis_expr = begin
                if haskey(_BASIS_CACHE, $new_ty_esc)
                    _BASIS_CACHE[$new_ty_esc]
                else
                    x = (() -> $(esc(expr)))()
                    x = convert(Vector{$base_ty_esc}, x)
                    @assert allunique(x)
                    _BASIS_CACHE[$new_ty_esc] = map($new_ty_esc, x)
                end
            end
        end
    end
end
