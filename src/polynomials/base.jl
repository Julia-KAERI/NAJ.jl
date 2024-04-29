# define abstract type for structure for polynomial types

abstract type AbstractPolynomial end
abstract type AbstractBasisPolynomial{T, X} end
abstract type PolynomialBase end

# Foundational functions for PolynomialBase

function (p::PolynomialBase)(x)
    return evalpoly(x, p.coeffs)
end

base_type(t::AbstractBasisPolynomial{T, X}) where {T, X} = X

# function (p::AbstractBasisPolynomial)(x)
#     bf = base_type(p) 
#     r = ([p.coeffs[i]* ((bf)(i-1))(x) for i in 1:length(p.coeffs)]) |> sum
#     return r
# end    

# Base.length(p::AbstractBasisPolynomial) = length(p.coeffs)

# order(p::AbstractBasisPolynomial) = length(p)-1
# degree(p::AbstractBasisPolynomial) = order(p)


# Foundational functions for AbastractBasisPolynomial

function (p::AbstractBasisPolynomial)(x)
    bf = base_type(p) 
    r = ([p.coeffs[i]* ((bf)(i-1))(x) for i in 1:length(p.coeffs)]) |> sum
    return r
end    

Base.length(p::AbstractBasisPolynomial) = length(p.coeffs)

order(p::AbstractBasisPolynomial) = length(p)-1
degree(p::AbstractBasisPolynomial) = order(p)
