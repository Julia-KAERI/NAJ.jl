function generate_legendre_polynomial_coefficients(n::Integer = 20)
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
        return coeffs
    end
end

function generate_chevyshev_polynomial_coefficients(n::Integer = 20)
    @assert n ≥ 0
    
    coeffs = [[1//1, ], [0//1, 1//1]]
    m = n+1
    
    if m < 2 
        return convert.(T, coeffs[m])
    else 
        for k in 3:m
            a1 = vcat([0], coeffs[k-1].*2) 
            a2 = vcat(coeffs[k-2], [0, 0], ) 
            push!(coeffs, a1.-a2 )
        end
        return coeffs
    end
end

function generate_hermite_polynomial_coefficients(n::Integer = 0)
    @assert n ≥ 0
    coeffs = [[1//1, ], [0//1, 2//1]]
    m = n+1
    
    if m < 3
        return coeffs
    else 
        for k in 3:m
            a1 = vcat([0], coeffs[k-1].*2) 
            a2 = vcat( (2 * (k-2)) .* coeffs[k-2], [0, 0], ) 
            push!(coeffs, a1.-a2 )
        end
        return coeffs
    end
end

function generate_laguerre_polynomial_coefficients(n::Integer = 0)
    @assert n ≥ 0
    coeffs = [[1//1, ], [1//1, -1//1]]
    m = n+1
    
    if m < 3
        return coeffs
    else 
        for k in 3:m
            a1 = vcat(coeffs[k-1].*(2*k-3), [0,])
            a2 = vcat([0], coeffs[k-1])
            a3 = vcat( (k-2) .* coeffs[k-2], [0, 0], ) 
            push!(coeffs, (a1 .- a2 .- a3).//(k-1) )
        end
        return coeffs
    end
end