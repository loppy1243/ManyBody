_s(s) = getfield(s, :_state)

struct Sub{B<:AbstractBasis, T} <: AbstractBasis
    _state::B

    function Sub{B, T}(b::B) where {B<:AbstractBasis, T}
        @assert hasmethod(Base.in, Tuple{B, Type{Sub{B, T}}})
        @assert b in Sub{B, T}

        new(b)
    end
end
Base.getproperty(s::Sub, prop::Symbol) = Base.getproperty(_s(s), prop)

superbasis(s::Type{<:Sub{B}}) where B<:AbstractBasis = B
superbasis(s::Sub) = superbasis(typeof(s))
supbasis(s::Type{<:Sub}) = superbasis(s)
supbasis(s::Type{<:Sub{<:Sub}}) = supbasis(superbasis(B))
supbasis(s::Sub) = supbasis(typeof(s))

superindex(s::Sub) = index(_s(s))
supindex(s::Sub) = superindex(s)
supindex(s::Sub{<:Sub}) = supindex(superindex(s))

convert(::Type{B}, s::Sub{B}) where B<:TensorBasis = _s(s)
convert(B::Type{<:TensorBasis}, s::Sub) = convert(B, _s(s))
B(s::Sub) where B<:TensorBasis = convert(B, s)

Base.:(==)(x::Sub, y::Sub) = _s(x) == _s(y)
Base.:(==)(x::TensorBasis, y::Sub) = x == _s(y)
Base.:(==)(x::Sub, y::TensorBasis) = _s(x) == y

Base.in(x, SB::Type{<:Sub}) = convert(superbasis(SB), x) in SB

let IXMAPS=Dict{Type{Sub}, Union{Vector{Int}, Vector{CartesianIndex}}}
    global subindexmap
    function subindexmap(SB::Sub)
        if !haskey(IXMAPS, SB)
            IXMAPS[SB] = [index(b) for b in superbasis(SB) if b in SB]
        end

        IXMAPS[SB]
    end
end

dim(SB::Type{<:Sub}) = length(subindexmap(SB))
index(sb::Sub) = findfirst(==(_s(sb)), subindexmap(typeof(sb)))
indexbasis(SB::Type{<:Sub}, i) = subindexmap(SB)[i]

macro defSub(f, ty_expr::Expr)
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

    in_expr(b) = add_where(:(Base.in($b::$base_ty_esc, ::Type{$new_ty_esc})))

    clos_param(s::Symbol) = s
    clos_param(s::Expr) = (@assert(s.head === :tuple && length(s.args) == 1); s.args[1])

    quote
        struct $subbasis_ty_esc end
        global const $new_ty_esc = Sub{$base_ty_esc, $subbasis_ty_esc}
        $(in_expr(clos_param(f.args[1]))) = $(esc(f.args[2]))
    end
end
