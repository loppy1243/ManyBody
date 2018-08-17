export normalize, rank

# Slater{B} == Product{N, NTuple{N, B}} where N for 0 <= N <= dim(B)

#abstract type AbstractProduct <: AbstractBasis
#
#struct BitProduct{N, BS<:NTuple{N, AbstractBasis}} <: AbstractProduct
#    bits::BitArray{N}
#end

struct Neg{B<:AbstractBasis} <: AbstractBasis
    state::B
end
const MaybeNeg{B<:AbstractBasis} = Union{B, Neg{B}}
#Neg(state::AbstractBasis) = Neg{typeof(state)}(state)
Base.:+(s::AbstractState) = s
Base.:-(b::AbstractBasis) = Neg(b)

innertype(::Type{Neg{B}}) where B<:AbstractBasis = B
inner(b::Neg) = b.state

index(b::Neg) = index(inner(b))
indexbasis(B::Type{<:Neg}, i::Int) = Neg(innertype(B)[i])
dim(B::Type{<:Neg}) = dim(innertype(B))

struct Product{N, BS<:NTuple{N, AbstractBasis}} <: AbstractBasis
    states::BS
end
Product(args::Vararg{AbstractBasis, N}) where N =
    Product{N, typeof(args)}(args)

rank(::Type{<:AbstractBasis}) = 1
rank(::Type{<:Product{N}}) where N = N

innertype(::Type{Product{1, Tuple{B}}}) where B<:AbstractBasis = B
inner(b::Product{1}) = b.states[1]

innertypes(::Type{Product{N, BS}}) where {N, BS<:NTuple{N, AbstractBasis}} =
    Tuple(BS.parameters)

## TODO: TESTME
@generated function normalize(b::Product)
    c = 0
    xs = map(innertypes(b)) do B
        if B <: Neg
            c += 1
            (:(inner(b.states[$i])), innertype(B))
        else
            (:(b.states[$i]), B)
        end
    end
    exprs, tys = zip(xs...)
    prod_expr = :(Product{$(rank(b)), Tuple{$(tys...)}}(($(exprs...),)))

    if iseven(c)
        prod_expr
    else
        :(Neg($prod_expr))
    end
end

## Note that we consider the products X*Y and Y*X as distinct (yet of course isomorphic).
Base.:(==)(a::B, b::B) where B<:Product = a.states == b.states

@generated function Base.promote_rule(P1::Type{<:Product{N}}, P2::Type{<:Product{N}}) where N
    ty_args = map(promote_type, innertypes(P1), innertypes(P2))

    :(Product{N, Tuple{$(ty_args...)}})
end
Base.convert(P1::Type{<:Product{N}}, p2::Product{N}) where N =
    P1(map(convert, innertypes(P1), p2.states))

Base.promote_rule(::Type{Product{1, Tuple{B}}}, ::Type{B}) where B<:AbstractBasis =
    Product{1, Tuple{B}}
Base.convert(::Type{Product{1, Tuple{B}}}, b::B) where B<:AbstractBasis = Product{1, Tuple{B}}(b)

## Should these exists? We lose a sign both ways.
Base.convert(::Type{Product{N, NTuple{N, B}} where N}, b::Slater{B}) where B<:AbstractBasis =
    Product(Tuple(B[i] for i = 1:dim(B) if b.bits[i]))
Base.convert(::Type{Slater{B}}, b::Product{N, NTuple{N, B}}) where N = Slater{B}(b.states)

innerindices(b::Product) = map(index, b.states)

index(b::Product) = index(b.states[1]) + sum(enumerate(b.states[2:end])) do I
    i, c = I
    (index(c)-1)*prod(dim.(typeof.(b.states[1:i])))
end

## Could @generate to eliminate func calls, unroll loop
function indexbasis(P::Type{<:Product}, i::Int)
    ixs = Vector{AbstractBasis}(undef, rank(P))

    i -= 1
    for (k, B) in enumerate(innertypes(P))
        j = i % dim(B) + 1
        ixs[k] = B[j]
        i = div(i - j + 1, dim(B))
    end

    P(Tuple(ixs))
end

innerdims(B::Type{<:Product}) = map(dim, innertypes(B))
dim(B::Type{<:Product}) = prod(dims(B))

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
    exprs = map(enumerate(bs)) do X
        i, b = X
        if b <: Product
            [:(bs[$i].states[$j]) for j = 1:rank(b)]
        elseif b <: Neg{<:Product}
            c += 1
            [:(inner(bs[$i]).states[$j]) for j = 1:rank(innertype(b))]
        elseif b <: Neg
            c += 1
            [:(inner(bs[$i]))]
        else
            [:(bs[$i])]
        end
    end |> x -> reduce(vcat, x)

    prod_expr = :(Product{$(length(bs)), $bs}(($(exprs...),)))
    if iseven(c)
        prod_expr
    else
        :(Neg($prod_expr))
    end
end
