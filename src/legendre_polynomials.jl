function legendre_polynomial_coefficients(T::Type=Float64, n::Integer = 0)
    @assert n ≥ 0
    
    coeffs = [[1//1, ], [0//1, 1//1]]
    m = n+1
    
    if m < 2 
        return coeffs[m]
    else 
        for k in 3:m
            a1 = vcat([0], coeffs[k-1]) .* ((2k-3)//(k-1))
            a2 = vcat(coeffs[k-2], [0, 0], ) .*((k-2)//(k-1))
            push!(coeffs, a1.-a2 )
        end
        return convert.(T, coeffs[end])
    end
end



struct LegendrePolynomial{T}<:AbstractSimplePolynomial
    order::Integer
    coeffs::Vector{T}

    function LegendrePolynomial(order::Integer = 0) 
        @assert order ≥ 0
        return new{Float64}(order, legendre_polynomial_coefficients(Float64, order))
    end

    function LegendrePolynomial{T}(order::Integer = 0) where T<:Real
        @assert order ≥ 0
        return new{T}(order, legendre_polynomial_coefficients(T, order))
    end
end

# order(p::LegendrePolynomial) = p.order
# degree(p::LegendrePolynomial) = order(p)

function (p::LegendrePolynomial)(x::T) where T <: Number
    return evalpoly(x, p.coeffs)
end


function Base.show(io::IO, p::LegendrePolynomial{T}) where T<:Real
    println(io, "LegendrePolynomial{$T}($(p.order))")
end

function toSimplePolynomial(p::LegendrePolynomial)
    return SimplePolynomial(p.coeffs)
end


