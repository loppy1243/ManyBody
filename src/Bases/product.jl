# Slater{B} == Product{NTuple{N, B}, Inf} where N for 0 <= N <= dim(B)

## Really should be Products of AbstractBasis
struct Product{M, BS<:NTuple{<:Any, TensorBasis}} <: TensorBasis{M}
    _states::BS

    function Product{M, BS}(states::BS) where {M, BS<:NTuple{<:Any, AbstractBasis}}
        @assert M == sum(rank, states)
        new(states)
    end
end
Product(args::Tuple) = Product{sum(rank, args), typeof(args)}(args)
Product(args::TensorBasis...) = Product(args)
subbases(::Type{<:Product{<:Any, BS}}) where BS<:NTuple{<:Any, TensorBasis} = BS.types

## Note that we consider the products X*Y and Y*X as distinct (yet of course isomorphic).
Base.:(==)(a::B, b::B) where B<:Product = a._states == b._states

index(b::Product) = CartesianIndex(map(index, b._states))
function indexbasis(B::Type{<:Product{M}}, ixs::Vararg{Int, M}) where M
    SBs = subbases(B)
    ret = Vector{TensorBasis}(undef, length(SBs))
    for (i, SB) in enumerate(SBs)
        R = rank(SB)
        ret[i] = indexbasis(SB, ixs[(i-1)*R+1:i*R]...)
    end

    B(Tuple(ret))
end
