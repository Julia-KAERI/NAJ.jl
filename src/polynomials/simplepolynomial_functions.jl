"""
    monic(p::P) where {P<:SimplePolynomial}

return monic polynomial of which highest coefficient is 1
"""
function monic(p::SimplePolynomial)
    return p/p.coeffs[end]
end

"""
    differentiate(p::SimplePolynomial)

return derivatives of polynomial p
"""
function differentiate(p::SimplePolynomial)
    if length(p) < 2 
        return SimplePolynomial([one(eltype(p.coeffs)), ])
    else
        coeffs = p.coeffs[2:end] .* (1:(length(p)-1))
        return SimplePolynomial(coeffs)
    end
end


"""
    integrate(p::P, a, b)

Integrate polynomial p. If both a and b are numbers, return definite integration from a to b. If a and b
are nothing, returns indefinite integral with constant 0. If only one of a and b are nothing, returns
indefinite integral with constant of non-nothing number.
"""
function integrate(p::SimplePolynomial, a::Union{Nothing, Number}=nothing, b::Union{Nothing, Number}=nothing) 

    if eltype(p.coeffs) <: Integer
        coeffs = zeros(Float64, length(p)+1)
    else 
        coeffs = zeros(eltype(p.coeffs), length(p)+1)
    end
    
    for i in 1:length(p.coeffs)
        coeffs[i+1] = p.coeffs[i]/(i)
    end
    
    
    if a === nothing && b === nothing # 상수항이 0 인 부정적분
        coeffs[1] = zero(eltype(coeffs))
        return SimplePolynomial(coeffs)
    elseif a === nothing || b === nothing # 상수항이 a 혹은 b 로 주어진 부정적분
        coeffs[1] = a
        return SimplePolynomial(coeffs)
    else # a 에서 b 구간 까지의 정적분
        return evalpoly(b, coeffs) - evalpoly(a, coeffs)
    end
end


"""
    polynomial_from_roots(xp::Vector{T}) where T<:Number

return monic polynomial having roots xp[1],..., xp[end].
"""
function polynomial_from_roots(xp::AbstractVector{T}) where T<:Number 
    return prod([SimplePolynomial([-x0, 1]) for x0 in xp])
end

"""
    valdermond_polynomial(xp, yp)

return Valdermond Polynomial. 
"""
function valdermond_polynomial(
    xp::AbstractVector{T1}, 
    yp::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    
    N = length(xp)
    @assert length(xp) == length(yp)
    V = [x^(j-1) for x in xp, j in 1:length(xp)]
    return SimplePolynomial(V\yp)
end

"""
    valdermond_polynomial(xp, yp)

return Lagrange Polynomial. 
"""
function lagrange_polynomial(
    xp::AbstractVector{T1}, 
    yp::AbstractVector{T2}) where {T1<:Number, T2<:Number}

    N = length(xp)
    @assert length(xp) == length(yp)
    
    r = SimplePolynomial([zero(T2), ])
    for i in 1:N
        coef = yp[i]
        rt = one(T2)
        for j in 1:N
            if i ≠ j
                @inbounds coef = coef/(xp[i]-xp[j])
                @inbounds rt = rt*SimplePolynomial([-xp[j], 1.0])
            end
        end
        r += rt*coef
    end
    return r
end

"""
    valdermond_polynomial(xp, yp)

return Newton Polynomial.
"""

function newton_polynomial(
    xp::AbstractVector{T1}, 
    yp::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    n = length(xp)    
    @assert n == length(yp)
    T = promote_type(T1, T2)
    N = LowerTriangular(ones(T, n, n))
    for j in 2:n, i in j:n
        @inbounds N[i, j] = N[i, j-1]*(xp[i] - xp[j-1]) 
    end
    a = N\yp
    r = SimplePolynomial([a[1], ])
    for i in 2:(n)
        @inbounds r += a[i] * polynomial_from_roots(xp[1:i-1])
    end
    return r

end

"""
    least_square_poly(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, order::Integer)

least square approximation for x and y with order n polynomial. Return SimplePolynomial object.
"""
function least_square_poly(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, order::Integer)
    @assert length(x) == length(y)
    @assert order ≥ 2 && order < length(x)-1
    X = [(xi)^i for xi in x, i in 0:order]
    a = inv(X'*X) *X' * y
    return SimplePolynomial(a)
end