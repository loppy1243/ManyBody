export SubBasis, subindices, subindexes
using Combinatorics: combinations
using Loppy.Util: cartesian_pow
using ..SpinMod

struct SubBasis{B<:Basis, T} <: Basis
    state::B
end

Base.convert(::Type{B}, x::SubBasis{B}) where B = x.state
Base.convert(::Type{C}, x::SubBasis{B}) where {C, B<:C} = x.state
Base.convert(::Type{SB}, x::B) where SB<:SubBasis{B} = SB(x)

Base.:(==)(x::SB, y::SB) where SB<:SubBasis = x.state == y.state
Base.in(p::SP, s::SubBasis{PartHole{R}}) where {SP, R<:RefState{SP}} = p in s.state

dim(::Type{SB}) where SB<:SubBasis = length(subindices(SB))
index(x::SubBasis) = findfirst(basis(typeof(x)), x)
indexbasis(::Type{SB}, i) where {B, SB<:SubBasis{B}} = basis(SB)[i]

subindices(::Type{B}) where B<:Basis = indices(B)
subindexes(::Type{SB}) where SB<:SubBasis = map(Index{SB}, subindices(SB))

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
        const $new_ty_esc = SubBasis{$base_ty_esc, $subbasis_ty_esc}

        @generated $basis_expr = begin
            x = (() -> $(esc(expr)))()
            x = convert(Vector{$base_ty_esc}, x)
            @assert allunique(x)
            map($new_ty_esc, x)
        end
        @generated $(subindices_expr(:T)) = map(index, basis(T.parameters[1]))
    end
end

@defSubBasis PHPaired{F, L} <: PartHole{Fermi{F, Pairing{L}}} begin
    R = Fermi{F, Pairing{L}}
    SP = Pairing{L}

    ret = [PartHole{R}()]
    for ph = 1:L-F
        for lhs in combinations(1:F, ph), lps in combinations(F+1:L, ph)
            pairs = map(lhs, lps) do lh, lp
                [(SP(lp, SPINUP),   SP(lh, SPINUP)  ),
                 (SP(lp, SPINDOWN), SP(lh, SPINDOWN)),
                 (SP(lp, SPINUP),   SP(lh, SPINDOWN)),
                 (SP(lp, SPINDOWN), SP(lh, SPINUP)  )]
            end |> x -> reduce(vcat, x)

            push!(ret, PartHole{R}(pairs))
        end
    end

    ret
end
