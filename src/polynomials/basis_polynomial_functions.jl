
function Base.zero(a::P) where P<:AbstractBasisPolynomial
    return P([zero(eltype(a.coeffs)), ])
end

function Base.one(a::P) where P<:AbstractBasisPolynomial
    return P([one(eltype(a.coeffs)), ])
end

function Base.:-(b::P) where {P<: AbstractBasisPolynomial}
    coeffs = -b.coeffs
    return P(coeffs)
end

function Base.:+(a::T, b::AbstractBasisPolynomial{P, X}) where {T <: Number, P <: Number, X} 
    rT = promote_type(T, P)
    coeffs = rT.(b.coeffs)
    coeffs[1] += a
    return (typeof(b))(coeffs)
end

function Base.:+(b::AbstractBasisPolynomial{P}, a::T) where {P <: Number, T <: Number} 
    return a+b
end

function Base.:+(a::SimplePolynomial{P1}, b::SimplePolynomial{P2}) where {P1 <: Number, P2 <: Number} 
    rT = promote_type(P1, P2)
    if length(b) > length(a)
        coeffs = zeros(rT, length(b))
        coeffs[1:length(a)] = a.coeffs[:]
        coeffs += b.coeffs
    else 
        coeffs = zeros(rT, length(a))
        coeffs[1:length(b)] = b.coeffs[:]
        coeffs += a.coeffs
    end
    return SimplePolynomial(coeffs)
end

function Base.:+(a::AbstractBasisPolynomial{P1}, b::AbstractBasisPolynomial{P2}) where {P1 <: Number, P2 <: Number} 
    if base_type(a) != base_type(b)
        return (SimplePolynomial(a) + SimplePolynomial(b))
    else 
        rT = promote_type(P1, P2)
        if length(b) > length(a)
            coeffs = zeros(rT, length(b))
            coeffs[1:length(a)] = a.coeffs[:]
            coeffs += b.coeffs
        else 
            coeffs = zeros(rT, length(a))
            coeffs[1:length(b)] = b.coeffs[:]
            coeffs += a.coeffs
        end
        return BasisPolynomial(coeffs, base_type(a))
    end
end


function Base.:-(a::AbstractBasisPolynomial{P1}, b::AbstractBasisPolynomial{P2}) where {P1 <: Number, P2 <: Number} 
    return a + (-b)
end

function Base.:-(b::AbstractBasisPolynomial{P}, a::T) where {P <: Number, T <: Number} 
    return b+(-a)
end

function Base.:-(a::T, b::AbstractBasisPolynomial{P}) where {T <: Number, P <: Number} 
    return a+(-b)
end

function Base.:*(a::T, b::AbstractBasisPolynomial{P}) where {T <: Number, P <: Number} 
    pT = promote_type(P, T)
    return BasisPolynomial(b.coeffs*a, base_type(b))
end

function Base.:*(b::AbstractBasisPolynomial{P}, a::T) where {P <: Number, T <: Number} 
    return a*b
end

function Base.:*(a::SimplePolynomial{P1}, b::SimplePolynomial{P2}) where {P1 <: Number, P2 <:Number} 
    rT = promote_type(P1, P2)
    ord1, ord2 = degree(a), degree(b)
    ord = ord1*ord2
    coef = zeros(rT, ord+2)
    
    for i in 0:ord1, j in 0:ord2
        @inbounds coef[i+j+1] += a.coeffs[i+1]*b.coeffs[j+1]
    end
    return SimplePolynomial(coef)
end

function Base.:*(a::AbstractBasisPolynomial{P1}, b::AbstractBasisPolynomial{P2}) where {P1 <: Number, P2 <:Number} 
    return SimplePolynomial(a)*SimplePolynomial(b)
end

function Base.:/(b::SimplePolynomial{P}, a::T) where {P <: Number, T <: Number} 
    return b*(1/a)
end

