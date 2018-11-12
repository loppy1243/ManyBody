# Slater{B} == Product{N, NTuple{N, B}} where N for 0 <= N <= dim(B)

#abstract type AbstractProduct <: ConcreteBasis
#
#struct BitProduct{N, BS<:NTuple{N, ConcreteBasis}} <: AbstractProduct
#    bits::BitArray{N}
#end

struct Product{BS<:NTuple{<:Any, ConcreteBasis}, M} <: TensorBasis{M}
    _states::BS

    function Product{M, BS}(states::BS) where {M, BS<:NTuple{<:Any, ConcreteBasis}}
        @assert M == sum(rank, states)
        new(states)
    end
end
Product(args::Tuple) = Product{sum(rank, args), typeof(args)}(args)
Product(args::TensorBasis...) = Product(args)

## Note that we consider the products X*Y and Y*X as distinct (yet of course isomorphic).
Base.:(==)(a::B, b::B) where B<:Product = a._states == b._states

index(b::Product) = CartesianIndex(map(index, b._states))
@generated function indexbasis(::Type{B}, ixs::Vararg{Int, M}) where
                              {BS<:NTuple{<:Any, TensorBasis}, M, B<:Product{BS, M}}
    ranks = vcat([0], map(rank, BS.types))
    ranges = [1+ranks[i-1]:ranks[i] for i=2:length(ranks)]

    ix_exprs = Expr(:tuple)
    for (T, range) in zip(BS.types, ranges)
        push!(ix_exprs.args, :($T[$(range...)]))
    end

    :(B($ix_exprs))
end
