struct Sub{B<:TensorBasis, M, T} <: TensorBasis{M}
    state::B
end
const MaybeSub{B<:ConcreteBasis} = Union{B, Sub{B}}

Base.:(==)(x::S, y::S) where S<:Sub = x.state == y.state
@commutes Base.:(==)(x::B, y::S) where B<:TensorBasis, S<:Sub{B} =
    x == y.state

fulldims

dim(B::Type{<:Sub}) = length(basis(B))
index(b::Sub) = findfirst(==(b), basis(typeof(b)))
indexbasis(B::Type{<:Sub}, ix::Union{Int, Base.CartesianIndex}) = basis(B)[ix]

innertype(::Type{<:Sub{B}}) where B<:ConcreteBasis = B
inner(s::Sub{<:ConcreteBasis}) = s.state

## Example
#@defSub Paired{P, L} <: Slater{Pairing{L}} do s
#    SP = Pairing{L}
#
#    P == count(findall(s.occ)) do I
#        s.occ[flipspin(SP[I])]
#    end
#end

### Update line
##############################################################################################

macro defSub(ty_expr::Expr, expr)
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
    dict_key_ty = add_where(:(Type{$new_ty_esc}))
    dict_val_ty = add_where(:(Vector{$new_ty_esc}))

    quote
        struct $subbasis_ty_esc end
        const $new_ty_esc = Sub{$base_ty_esc, $subbasis_ty_esc}

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
