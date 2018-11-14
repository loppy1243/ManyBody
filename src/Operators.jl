@reexport module Operators
export tabulate, tabulate!, normord, applyop, applyop!
using ..ManyBody

### Tabulation
##############################################################################################
function tabulate!(f, A::AbstractArray, BL::Type{<:AbstractBasis}, BR::Type{<:AbstractBasis}
                  ;leftindices=nothing, rightindices=nothing)
    Il = leftindices === nothing ? tabulate_lindices(typeof(A), BL) : leftindices
    Ir = rightindices === nothing ? tabulate_rindices(typeof(A), BR) : rightindices
    for (bl, br) in Iterators.product(BL, BR)
        _tabulate!_kernel(A, Il, bl, Ir, br, f)
    end

    A
end
_tabulate!_kernel(A, Il, bl, Ir, br, f) = A[Il[bl], Ir[br]] = f

tabulate_lindices(A, B) = tabulate_indices(A, B)
tabulate_rindices(A, B) = tabulate_indices(A, B)
tabulate_indices(::Type{<:AbstractArray}, B) = CartesianIndices(B)
tabulate_indices(::Type{<:AbstractMatrix}, B) = LinearIndices(B)

tabulate(f, T, B) = tabulate(f, T, B, B)
tabulate(f, T::Type, BL, BR) = tabulate(f, Array{T}, BL, BR)
@generated function tabulate(f, ::Type{A}, ::Type{BL}, ::Type{BR}) where
                 {A<:AbstractArray, BL<:AbstractBasis, BR<:AbstractBasis}
    lindices = tabulate_lindices(A, BL)
    rindices = tabulate_rindices(A, BR)
    dims = (size(lindices)..., size(rindices)...)

    :(tabulate!(f, similar(A, $dims), BL, BR; leftindices=lindices, rightindices=rindices))
end

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

const AnyRLOp{SP<:AbstractBasis} = Union{RaiseLowerOp{SP}, RaiseOp{SP}, LowerOp{SP}}

(op::RaiseOp)(args...) = applyop(op, args...)
(op::LowerOp)(args...) = applyop(op, args...)
(op::RaiseLowerOp)(args...) = applyop(op, args...)

function applyop!(op::RaiseOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis
    sgn = create!(b, op.state)
    if sgn > 0
        (1, b)
    elseif sgn == 0
        (0, b)
    else
        (-1, b)
    end
end
function applyop!(op::LowerOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis
    sgn = annihil!(b, op.state)
    if sgn > 0
        (1, b)
    elseif sgn == 0
        (0, b)
    else
        (-1, b)
    end
end
function applyop(op::RLOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis
    ret = deepcopy(b)
    applyop!(op, ret)
end

applyop!(op::RaiseLowerOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis =
    foldr(applyop!, op.rlops, init=b)
function applyop(op::RaiseLowerOp{SP}, b::Bases.Slater{SP}) where SP<:AbstractBasis
    ret = deepcopy(b)
    applyop!(op, ret)
end

struct Raised{T}; val::T end
↑(x) = Raised(x)

macro A(p_exprs...)
    raised_exprs = map(p_exprs) do p_expr
        if p_expr isa Expr && p_expr.head == Symbol("'")
            :(Raised($(p_expr.args[1])))
        else
            p_expr
        end
    end

    :(_A($(raised_exprs...)))
end

_A(p::AbstractBasis) = LowerOp(p)
_A(p::Raised{<:AbstractBasis}) = RaiseOp(p.val)
@generated function _A(sps::Vararg{Union{B, Raised{B}}}) where B<:AbstractBasis
    args = map(enumerate(sps)) do x
        i, U = x
        if U <: Raised
            :(RaiseOp(sps[$i].val))
        else
            :(LowerOp(sps[$i]))
        end
    end

    :(RaiseLowerOp{B}(RLOp{B}[$(args...)]))
end

normord(a::RaiseLowerOp) = normord(RefStates.Vacuum(spbasis(a)), a)
function normord(R::RefState{B}, a::RaiseLowerOp{B}) where B<:AbstractBasis
    function comp(a, b)
        if a[2] && b[2]
            index(a[1]) >= index(b[1])
        elseif !a[2] && !b[2]
            index(a[1]) <= index(b[1])
        else
            a[2]
        end
    end

    xs = map(a.rlops) do b
        (b.state, if b isa RaiseOp
            ~isocc(R, b.state)
        else
            isocc(R, b.state)
        end)
    end

    p = sortperm(xs, lt=comp)
    sgn = levicivita(p)
    ret = RaiseLowerOp{B}(a.rlops[p])

    iszero(sgn) && error("IMPOSSIBLE")

    (sgn, ret)
end

Base.show(io::IO, x::RaiseOp) = print(io, "A($(x.state)')")
Base.show(io::IO, x::LowerOp) = print(io, "A($(x.state) )")
function Base.show(io::IO, x::RaiseLowerOp)
    f(a::RaiseOp) = "$(a.state)'"
    f(a::LowerOp) = "$(a.state) "

    print(io, "A(", f(x.rlops[1]))
    for op in x.rlops[2:end]
        print(io, ", $(f(op))")
    end
    print(io, ")")
end

end # module Operators
