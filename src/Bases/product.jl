# Slater{B} == Product{N, NTuple{N, B}} where N for 0 <= N <= dim(B)

#abstract type AbstractProduct <: AbstractBasis
#
#struct BitProduct{N, BS<:NTuple{N, AbstractBasis}} <: AbstractProduct
#    bits::BitArray{N}
#end

struct Neg{B<:AbstractBasis} <: AbstractBasis
    state::B
end
Neg(state::AbstractBasis) = Neg{typeof(state)}(state)
Base.:+(b::AbstractBasis) = b
Base.:-(b::AbstractBasis) = Neg(b)

innertype(::Type{Neg{B}}) where B<:AbstractBasis = B
inner(b::Neg) = b.state

index(b::Neg) = index(inner(b))
indexbasis(B::Type{<:Neg}, i::Int) = Neg(innertype(B)[i])
dim(B::Type{<:Neg}) = dim(innertype(B))

struct Product{N, BS<:NTuple{N, AbstractBasis}} <: AbstractBasis
    states::BS
end
const CF64Product{N, BS<:NTuple{N, AbstractBasis}} = Product{N, BS, ComplexF64}
Product(args::Vararg{AbstractBasis, N}) where N =
    Product{N, typeof(args)}(args)

Base.convert(::Type{Product{1, Tuple{B}}}, b::B) = Product{1, Tuple{B}}(b)
Base.promote(::Type{Product{1, Tuple{B}}}, ::Type{B}) where B<:AbstractBasis =
    Product{1, Tuple{B}}

index(b::Product) = prod(enumerate(b.states)) do I
    i, c = I
    index(c)*prod(dim.(typeof.(b.states[1:i-1])))
end

_states(a::AbstractBasis) = (a,)
_states(a::Product) = a.states
function Base.:*(bs::AbstractBasis...)
    ret = AbstractBasis[]
    for b in bs
        for s in _states(b)
            push!(ret, s)
        end
    end

    Product(Tuple(ret))
end
