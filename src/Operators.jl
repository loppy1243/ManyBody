@reexport module OperatorsMod
export Operator, tabulate

using ..States

struct Operator{B, Basis<:State, Op}
    op::Op
end

(op::Operator{B, Basis, <:Function})(sps::Vararg{Basis, 2B}) where {B, Basis} = op.op(sps...)
(op::Operator{B, Basis, <:AbstractArray{T, 2B}})(sps::Vararg{Basis, 2B}) where {B, Basis, T} =
    op.op[map(pnum, sps)...]

tabulate(op::Operator) = tabulate(Float64, op)
tabulate(op::Operator{B, Basis, <:AbstractArray}) where {B, Basis} = op
function tabulate(T::Type, op::Operator{B, Basis, Op}) where {B, Basis, Op}
    states = iter(Basis) |> collect
    statenums = map(pnum, states)
    numstates = length(states)

    mat = Array{T, 2B}(fill(numstates, 2B)...)
    for ixs in Iterators.product(fill(1:numstates, 2B)...)
        ixs = collect(ixs)
        mat[statenums[ixs]...] = op(states[ixs]...)
    end

    Operator{B, Basis, typeof(mat), T}(mat)
end

end # module OperatorsMod
