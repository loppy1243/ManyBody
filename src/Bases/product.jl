# Slater{B} == Product{N, NTuple{N, B}} where N for 0 <= N <= dim(B)

#abstract type AbstractProduct <: AbstractBasis
#
#struct BitProduct{N, BS<:NTuple{N, AbstractBasis}} <: AbstractProduct
#    bits::BitArray{N}
#end

struct Neg{B<:AbstractBasis} <: AbstractBasis
    state::B
end
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

## Note that we consider the products X*Y and Y*X as distinct (yet of course isomorphic).
Base.:(==)(a::B, b::B) where B<:Product = a.states == b.states

Base.convert(::Type{Product{1, Tuple{B}}}, b::B) where B<:AbstractBasis = Product{1, Tuple{B}}(b)
Base.promote_rule(::Type{Product{1, Tuple{B}}}, ::Type{B}) where B<:AbstractBasis =
    Product{1, Tuple{B}}

index(b::Product) = index(b.states[1]) + sum(enumerate(b.states[2:end])) do I
    i, c = I
    (index(c)-1)*prod(dim.(typeof.(b.states[1:i])))
end

## Could @generate to eliminate func calls, unroll loop
function indexbasis(::Type{Product{N, BS}}, i::Int) where {N, BS<:NTuple{N, AbstractBasis}}
    ixs = Vector{AbstractBasis}(undef, N)

    i -= 1
    for (k, B) in enumerate(BS.parameters)
        j = i % dim(B) + 1
        ixs[k] = B[j]
        i = div(i - j + 1, dim(B))
    end

    Product{N, BS}(Tuple(ixs))
end

dim(::Type{Product{N, BS}}) where {N, BS<:NTuple{N, AbstractBasis}} =
    prod(dim, BS.parameters)
rank(::Type{<:Product{N}}) where N = N

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
    exprs = map(enumerate(bs)) do X
        i, b = X
        if b <: Product
            [:(bs[$i].states[$j]) for j = 1:rank(b)]
        else
            [:(bs[$i])]
        end
    end |> x -> reduce(vcat, x)

    :(Product{$(length(bs)), $bs}(($(exprs...),)))
end
