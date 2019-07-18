@reexport module Operators
export normord, applyop, applyop!
using ..ManyBody
using Combinatorics: levicivita
using Base.Cartesian

### Raising and Lowering Operators
##############################################################################################
struct RaiseOp{SP<:AbstractBasis}
    state::SP
end
struct LowerOp{SP<:AbstractBasis}
    state::SP
end
const RLOp{SP<:AbstractBasis} = Union{RaiseOp{SP}, LowerOp{SP}}

struct RaiseLowerOp{SP<:AbstractBasis}
    rlops::Vector{RLOp{SP}}

    RaiseLowerOp{SP}(rlops::Vector{RLOp{SP}}) where SP<:AbstractBasis = new(rlops)
end
RaiseLowerOp{SP}(itr) where SP<:AbstractBasis = RaiseLowerOp(collect(RLOp{SP}, itr))
RaiseLowerOp(rlops::Vector{RLOp{SP}}) where SP<:AbstractBasis = RaiseLowerOp{SP}(rlops)

Base.:(==)(op1::O, op2::O) where O<:RLOp = op1.state == op2.state
Base.:(==)(op1::RLOp, op2::RLOp) = false

const AnyRLOp{SP<:AbstractBasis} = Union{RaiseLowerOp{SP}, RaiseOp{SP}, LowerOp{SP}}
spbasis(::Type{<:AnyRLOp{SP}}) where SP<:AbstractBasis = SP
spbasis(a::AnyRLOp) = spbasis(typeof(a))

@generated (op::RaiseOp)(args::Vararg{Any, N}) where N =
    :(@ncall($N, applyop, op, i -> args[i]))
@generated (op::LowerOp)(args::Vararg{Any, N}) where N =
    :(@ncall($N, applyop, op, i -> args[i]))
@generated (op::RaiseLowerOp)(args::Vararg{Any, N}) where N =
    :(@ncall($N, applyop, op, i -> args[i]))

_apply_rlop!(op::RaiseOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis = create!(b, op.state)
_apply_rlop!(op::LowerOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis = annihil!(b, op.state)
_apply_rlop(op::RaiseOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis = create(b, op.state)
_apply_rlop(op::LowerOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis = annihil(b, op.state)

for (opapp, app) in zip((:applyop!, :applyop), (:_apply_rlop!, :_apply_rlop))
    @eval function $opapp(op::RLOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis
        sgn = $app(op, b)[1]
        (sgn, b)
    end
end

function applyop!(op::RaiseLowerOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis
    sgn = 1
    for rlop in op.rlops
        sgn *= applyop!(rlop, b)[1]
        if iszero(sgn)
            return (sgn, b)
        end
    end

    (sgn, b)
end
applyop(op::RaiseLowerOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis =
    applyop!(op, copy(b))

### Matrix elements
##############################################################################################
applyop(op::AnyRLOp, f::DualBasis, b::AbstractBasis) = applyop(op, f', b)

@generated applyop(op::RLOp{SP}, b1::Bases.Slater{SP}, b2::Bases.Slater{SP}) where
                   SP<:AbstractBasis = quote
    sgn = acsign(b2, op.state)
    sgn*$(op<:RaiseOp ? :(~b2.bits[+op.state]) : :(b2.bits[+op.state]))
end

function applyop(op::RaiseLowerOp{SP}, b1::Bases.Slater{SP}, b2::Bases.Slater{SP}) where
                SP<:AbstractBasis
    true_bits = BitSet()
    false_bits = BitSet()

    ixs = Int[-x.state for x in op.rlops]

    sgn = 1
    for (k, (i, rlop)) in zip(ixs, op.rlops) |> enumerate
        sgn′, bit = _apply_rlop!(rlop, b2)
        if !(i in true_bits || i in false_bits)
            bit ? union!(true_bits,  i) : union!(false_bits, i)
        end
        sgn *= sgn′
#        # If we break here, then we have to account for that in true_bits and false_bits
#        if iszero(sgn)
#            resize!(ixs, k)
#            break
#        end
    end

    ret = iszero(sgn) ? sgn : sgn*overlap(b1, b2)

    for (i, rlop) in zip(ixs, op.rlops)
        # Note that i *must* be in one or the other
        b2.bits[i] = (i in true_bits) || ~(i in false_bits)
    end

    ret
end

### @A syntax
##############################################################################################
struct Raised{T}; val::T end

macro A(p_exprs...)
    raised_exprs = map(p_exprs) do p_expr
        if p_expr isa Expr && p_expr.head == Symbol("'")
            :(Raised($(esc(p_expr.args[1]))))
        else
            esc(p_expr)
        end
    end

    :(_A($(raised_exprs...)))
end

_A(p::AbstractBasis) = LowerOp(p)
_A(p::Raised{<:AbstractBasis}) = RaiseOp(p.val)
_A(sps::Vararg{Union{B, Raised{B}}}) where B<:AbstractBasis =
    RaiseLowerOp{B}(reverse(map(_A, sps)))

##############################################################################################
##############################################################################################
normord(a::RaiseLowerOp) = normord(RefStates.Vacuum{spbasis(a)}(), a)
function normord(R::RefState{B}, a::RaiseLowerOp{B}) where B<:AbstractBasis
    function comp(a, b)
        if a isa RaiseOp && b isa LowerOp
            isunocc(R, a.state) || isunocc(R, b.state)
        elseif a isa LowerOp && b isa RaiseOp
            isocc(R, a.state) || isocc(R, b.state)
        elseif a isa RaiseOp && b isa RaiseOp
            if isocc(R, a.state) == isocc(R, b.state)
                index(a.state) <= index(b.state)
            else
                isunocc(R, a.state) || isocc(R, b.state)
            end
        else#if a isa LowerOp && b isa LowerOp
            if isocc(R, a.state) == isocc(R, b.state)
                index(a.state) > index(b.state)
            else
                isocc(R, a.state) || isunocc(R, b.state)
            end
        end
    end

    p = sortperm(a.rlops, lt=!comp)
    sgn = levicivita(p)
    ret = RaiseLowerOp{B}(a.rlops[p])

    @assert !iszero(sgn)

    (sgn, ret)
end

Base.show(io::IO, x::RaiseOp) = print(io, "A($(x.state)')")
Base.show(io::IO, x::LowerOp) = print(io, "A($(x.state) )")
function Base.show(io::IO, x::RaiseLowerOp)
    f(a::RaiseOp) = "$(a.state)'"
    f(a::LowerOp) = "$(a.state) "

    print(io, "A(", f(x.rlops[end]))
    for i in length(x.rlops)-1:-1:1
        print(io, ", $(f(x.rlops[i]))")
    end
    print(io, ")")
end

end # module Operators
