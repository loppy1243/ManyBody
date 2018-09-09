# Slater{B} == Product{N, NTuple{N, B}} where N for 0 <= N <= dim(B)

#abstract type AbstractProduct <: ConcreteBasis
#
#struct BitProduct{N, BS<:NTuple{N, ConcreteBasis}} <: AbstractProduct
#    bits::BitArray{N}
#end

struct Neg{B<:ConcreteBasis} <: LeafBasis
    state::B
end
const MaybeNeg{B<:ConcreteBasis} = Union{B, Neg{B}}
Base.:-(b::ConcreteBasis) = Neg(b)
Base.:-(b::Neg) = inner(b)

IndexType(B::Type{<:Neg}) = IndexType(innertype(B))

innertype(::Type{Neg{B}}) where B<:AbstractBasis = B
inner(b::Neg) = b.state

index(b::Neg) = index(inner(b))
indexbasis(B::Type{<:Neg}, ixs...) = Neg(innertype(B)[ixs...])
dim(B::Type{<:Neg}) = dim(innertype(B))
innerdims(B::Type{<:Neg}) = innerdims(inertype(B))

Base.convert(::Type{B}, b::Neg{B}) where B<:ConcreteBasis = inner(b)
Base.convert(B::Type{<:Neg}, b::ConcreteBasis) = B(b)

struct Product{N, BS<:NTuple{N, ConcreteBasis}} <: ConcreteBasis
    states::BS

    function Product{N, BS}(states::BS) where {N, BS<:NTuple{N, ConcreteBasis}}
        new(states)
    end
end
Product(args::Tuple) = Product{length(args), typeof(args)}(args)
Product(args::ConcreteBasis...) = Product(args)

IndexType(::Type{<:Product{N}}) where N = IndexTypes.Cartesian{N}()
rank(::Type{<:Product{N}}) where N = N

innertype(::Type{Product{1, Tuple{B}}}) where B<:AbstractBasis = B
inner(b::Product{1}) = b.states[1]

innertypes(::Type{Product{N, BS}}) where {N, BS<:NTuple{N, AbstractBasis}} =
    Tuple(BS.parameters)

index(b::Product) = Base.CartesianIndex(map(index, b.states))
indexbasis(P::Type{<:Product{N}}, ixs::Base.CartesianIndex{N}) where N = indexbasis(P, Tuple(ixs))
indexbasis(P::Type{<:Product{N}}, ixs::Vararg{Int, N}) where N = indexbasis(P, ixs)
indexbasis(P::Type{<:Product{N}}, ixs::NTuple{N, Int}) where N =
    P(map(indexbasis, innertypes(P), ixs))

innerdims(B::Type{<:Product}) = map(dim, innertypes(B))
dim(B::Type{<:Product}) = prod(innerdims(B))

## Note that we consider the products X*Y and Y*X as distinct (yet of course isomorphic).
Base.:(==)(a::B, b::B) where B<:Product = a.states == b.states

Base.getindex(b::Product, i::Int) = b.states[i]
Base.getindex(b::Product, ixs::AbstractArray) = map(i -> b.states[i], ixs)

Base.iterate(b::Product) = iterate(b.states)
Base.iterate(b::Product, st) = iterate(b.states, st)
Base.IteratorSize(B::Type{<:Product}) = Base.HasLength()
Base.length(b::Product) = rank(b)

#Base.convert(::Type{P}, p::P) where P<:Product = p
#Base.convert(::Type{P}, p::P) where P<:Product{1} = p

function Base.promote_rule(P1::Type{<:Product{N}}, P2::Type{<:Product{N}}) where N
    tys = map(promote_type, innertypes(P1), innertypes(P2))
    Product{N, Tuple{ty_args...}}
end
Base.convert(P1::Type{<:Product{N}}, p2::Product{N}) where N =
    P1(map(convert, innertypes(P1), p2.states))

Base.promote_rule(::Type{Product{1, Tuple{B}}}, ::Type{B}) where B<:ConcreteBasis =
    Product{1, Tuple{B}}
Base.convert(::Type{Product{1, Tuple{B}}}, b::B) where B<:ConcreteBasis = Product(b)
Base.convert(::Type{B}, b::Product{1, Tuple{B}}) where B<:ConcreteBasis = inner(b)

## Would an explicit tensor() function be better?

#_states(a::AbstractBasis) = (a,)
#_states(a::Product) = a.states
#function Base.:*(bs::AbstractBasis...)
#    ret = AbstractBasis[]
#    for b in bs
#        for s in _states(b)
#            push!(ret, s)
#        end
#    end
#
#    Product(Tuple(ret))
#end
@generated function Base.:*(bs::AbstractBasis...)
    c = 0
    tys, exprs = map(enumerate(bs)) do X
        i, b = X
        if b <: Product
            ([innertypes(b)...], [:(bs[$i].states[$j]) for j = 1:rank(b)])
        elseif b <: Neg{<:Product}
            c += 1
            ([innertypes(innertype(b))...],
             [:(inner(bs[$i]).states[$j]) for j = 1:rank(innertype(b))])
        elseif b <: Neg
            c += 1
            ([innertype(b)], [:(inner(bs[$i]))])
        elseif b <: LeafBasis
            ([innertype(b)], [:(inner(bs[$i]))])
        else
            ([b], [:(bs[$i])])
        end
    end |> x -> zip(x...)
    tys = reduce(vcat, tys)
    exprs = reduce(vcat, exprs)

    prod_expr = :(Product{$(length(bs)), Tuple{$(tys...)}}(($(exprs...),)))
    if iseven(c)
        prod_expr
    else
        :(Neg($prod_expr))
    end
end
