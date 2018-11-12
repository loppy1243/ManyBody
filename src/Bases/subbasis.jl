_s(s) = getfield(s, :_state)

struct Sub{B<:TensorBasis, M, T} <: TensorBasis{M}
    _state::B

    function Sub{B, M, T}(b::B) where {B<:TensorBasis, M, T}
        @assert isdefined(Bases.:subindexmap) #=
             =# && hasmethod(subindexmap, Tuple{Type{Sub{B, M, T}}})
        @assert findfirst(subindexmap(Sub{B, M, T}), index(b)) !== nothing

        new(b)
    end
end
Base.getproperty(s::Sub, prop::Symbol) = Base.getproperty(_s(s), prop)

superbasis(s::Sub{B}) where B<:TensorBasis = B
supbasis(s::Sub) = superbasis(s)
supbasis(s::Sub{<:Sub}) = supbasis(superbasis(B))

convert(::Type{B}, s::Sub{B}) where B<:TensorBasis = _s(s)
convert(B::Type{<:TensorBasis}, s::Sub) = convert(B, _s(s))
B(s::Sub) where B<:TensorBasis = convert(B, s)

Base.:(==)(x::Sub, y::Sub) = _s(x) == _s(y)
Base.:(==)(x::TensorBasis, y::Sub) = x == _s(y)
Base.:(==)(x::Sub, y::TensorBasis) = _s(x) == y

fulldims(B::Type{<:Sub}) = map(lastindex, axes(subindexmap(B)))
# Could store an inverse subindexmap as well.
index(b::Sub) = findfirst(==(index(_s(b))), subindexmap(typeof(b)))
indexbasis(SB::Type{<:Sub{<:TensorBasis, M}}, ixs::Vararg{Int, M}) =
    superbasis(SB)[subindexmap(SB)[CartesianIndex(ixs)]]

macro defSub(f, M, ty_expr::Expr)
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

    subindexmap_expr = add_where(:(Bases.subindexmap(::Type{$new_ty_esc})))

    quote
        let IXMAP = _makeixmap($f_esc, $base_ty_esc)
            struct $subbasis_ty_esc end
            global const $new_ty_esc = Sub{$base_ty_esc, ndims(IXMAP), $subbasis_ty_esc}

            $subindexmap_expr = IXMAP
        end
    end
end

function _makeixmap(f, B)
    N = rank(B)
    mask = map(f, eachindex(B))

    sizes = []
    for d = 1:N
        sz = count(any(mask, dims=[x for x=1:N if x != d]))
        !iszero(sz) && push!(sizes, zs)
    end

    ixmap = Array{CartesianIndices}(undef, sizes...)

    for (I, J) in zip(eachindex(ixmap), findall(mask))
        ixmap[I] = J
    end

    ixmap
end
