struct Sub{B<:AbstractBasis, T} <: AbstractBasis
    _state::B

    function Sub{B, T}(b::B) where {B<:AbstractBasis, T}
        isabstracttype(B) && @warn "Bases.Sub{B} defined with abstract B!"
        @assert hasmethod(Base.in, Tuple{B, _Arrayish{Sub{B, T}}})
        @assert b in elems(Sub{B, T})

        new(b)
    end
end
Base.getproperty(s::Sub, prop::Symbol) = Base.getproperty(superelem(s), prop)

superbasis(s::Type{<:AbstractBasis}) = s
superbasis(s::Type{<:Sub{B}}) where B<:AbstractBasis = B
superbasis(s::AbstractBasis) = superbasis(typeof(s))
supbasis(s::Type{<:AbstractBasis}) = s
supbasis(s::Type{<:Sub}) = superbasis(s)
supbasis(s::Type{<:Sub{<:Sub}}) = supbasis(superbasis(B))
supbasis(s::AbstractBasis) = supbasis(typeof(s))

superindex(s::AbstractBasis) = index(s)
superindex(s::Sub) = index(superelem(s))
supindex(s::AbstractBasis) = index(s)
supindex(s::Sub) = superindex(s)
supindex(s::Sub{<:Sub}) = supindex(superindex(s))

superelem(s::AbstractBasis) = s
superelem(s::Sub) = getfield(s, :_state)
supelem(s::AbstractBasis) = s
supelem(s::Sub) = superelem(s)
supelem(s::Sub{<:Sub}) = supelem(superelem(s))

Base.convert(::Type{B}, s::B) where B<:Sub = s
Base.convert(::Type{B}, s::Sub{B}) where B<:AbstractBasis = superelem(s)
Base.convert(B::Type{<:AbstractBasis}, s::Sub) = convert(B, superelem(s))
B(s::Sub) where B<:AbstractBasis = convert(B, s)

Base.:(==)(x::Sub, y::Sub) = supelem(x) == supelem(y)
Base.:(==)(x::AbstractBasis, y::Sub) = x == supelem(y)
Base.:(==)(x::Sub, y::AbstractBasis) = supelem(x) == y

Base.in(x, a::_Arrayish{<:Sub}) = convert(superbasis(eltype(a)), x) in elems(SB)

let IXMAPS=Dict{Any, Any}()
    global subindexmap
    function subindexmap(SB::Type{<:Sub})
        if !haskey(IXMAPS, SB)
            IXMAPS[SB] = [index(b) for b in elems(supbasis(SB)) if b in elems(SB)]
        end

        IXMAPS[SB]
    end
end

dim(SB::Type{<:Sub}) = length(subindexmap(SB))
index(sb::Sub) = findfirst(==(supindex(sb)), subindexmap(typeof(sb)))
indexbasis(SB::Type{<:Sub}, i::Int) = SB((+supbasis(SB))[subindexmap(SB)[i]])

macro defSub(ty_expr::Expr, body_expr::Expr)
    get_tys(s) = [s]
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
    @assert body_expr.head == :block #=
         =# && (body_expr.args[2] isa Symbol #=
             =# || body_expr.args[2] isa Expr && bod_expr.args[2].head == :tuple)

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

    in_expr(b) = add_where(:(Base.in($b::$base_ty_esc, ::_Arrayish{$new_ty_esc})))

    func_param_esc = esc(body_expr.args[2])
    body_esc = esc(Expr(:block, body_expr.args[3:end]...))

    quote
        struct $subbasis_ty_esc end
        global const $new_ty_esc = Sub{$base_ty_esc, $subbasis_ty_esc}
        $(in_expr(func_param_esc)) = $body_esc
    end
end

_depth(x::Type{<:Sub}, cnt) = _depth(superbasis(x), cnt+1)
_depth(x::Type{<:AbstractBasis}, cnt) = cnt
_depth(x::AbstractBasis) = _depth(typeof(x), 0)
Base.show(io::IO, x::Sub) = print(io, "Sub<$(_depth(x))>($(supelem(x)))")

function Base.show(io::IO, mime::MIME"text/plain", x::Sub)
    d = _depth(x)
    if d == 1
        println(io, typeof(x))
    else
        println(io, "Depth $(_depth(x)) $(typeof(x))")
    end
    show(io, mime, supelem(x))
end
