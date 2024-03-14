function (p::PolynomialBasis)(x)
    return evalpoly(x, p.coeffs)
end

basis_type(t::AbstractBasisPolynomial{T, X}) where {T, X} = X

function (p::AbstractBasisPolynomial)(x)
    bf = basis_type(p) 
    r = ([p.coeffs[i]* ((bf)(i-1))(x) for i in 1:length(p.coeffs)]) |> sum
    return r
end    

Base.length(p::AbstractBasisPolynomial) = length(p.coeffs)

order(p::AbstractBasisPolynomial) = length(p)-1
degree(p::AbstractBasisPolynomial) = order(p)

struct SimpleBasis{T}<:PolynomialBasis
    order::Integer
    coeffs::Vector{T} 

    function SimpleBasis(order::Integer)
        @assert order ≥ 0
        r = zeros(Float64, order+1)
        r[end]=one(Float64)
        return new{Float64}(order, r)
    end

    function SimpleBasis{T}(order::Integer) where T<:Real
        @assert order ≥ 0
        r = zeros(T, order+1)
        r[end]=one(T)
        return new{T}(order, r)
    end
end


"""
    legendre_polynomial_coefficients(T::Type=Float64, n::Integer = 0)
"""
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

struct LegendreBasis{T}<:PolynomialBasis
    order::Integer
    coeffs::Vector{T} 

    function LegendreBasis(order::Integer = 0) 
        @assert order ≥ 0
        return new{Float64}(order, legendre_polynomial_coefficients(Float64, order))
    end

    function LegendreBasis{T}(order::Integer = 0) where T<:Real
        @assert order ≥ 0
        return new{T}(order, legendre_polynomial_coefficients(T, order))
    end
end


function chevyshev_polynomial_coefficients(T::Type=Float64, n::Integer = 0)
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
        return convert.(T, coeffs[end])
    end
end

struct ChevyshevBase{T}<:PolynomialBasis
    order::Integer
    coeffs::Vector{T}

    function ChevyshevBase(order::Integer = 0) 
        @assert order ≥ 0
        return new{Float64}(order, chevyshev_polynomial_coefficients(Float64, order))
    end

    function ChevyshevBase{T}(order::Integer = 0) where T<:Real
        @assert order ≥ 0
        return new{T}(order, chevyshev_polynomial_coefficients(T, order))
    end
end

struct BasisPolynomial{T, X} <:AbstractBasisPolynomial{T, X}
    coeffs::Vector{T}
    # BaseType::X

    function BasisPolynomial(a::AbstractVector{P}, t::X) where {P <: Number, X<:Type}
        if length(a) == 0 
            return new{P, t}(zeros(T, 1))#, SimpleBasis{Float64})
        else 
            last_nz = findlast(!iszero, a)
            a_last = max(1, isnothing(last_nz) ? 0 : last_nz)
            return new{P, t}(a[1:a_last])#, SimpleBasis{Float64})
        end
    end

    function BasisPolynomial{T}(a::AbstractVector{P}, t::X) where {T <: Number, P<:Number, X<:Type}
        if length(a) == 0 
            return new{T, t}(zeros(T, 1))#, SimpleBasis{Float64})
        else 
            last_nz = findlast(!iszero, a)
            a_last = max(1, isnothing(last_nz) ? 0 : last_nz)
            return new{T, t}(convert.(T, a[1:a_last]))#, SimpleBasis{Float64})
        end
    end

end


const SimplePolynomial{T} = BasisPolynomial{T, SimpleBasis}
const LegendrePolynomial{T} = BasisPolynomial{T, LegendreBasis}
const ChevyshevPolynomial{T} = BasisPolynomial{T, ChevyshevBase}

function SimplePolynomial(a::AbstractVector{T}) where {T<:Real}
    return BasisPolynomial(a, SimpleBasis)
end

function LegendrePolynomial(a::AbstractVector{T}) where T<:Real
    return BasisPolynomial(a, LegendreBasis)
end

function ChevyshevPolynomial(a::AbstractVector{T}) where T<:Real
    return BasisPolynomial(a, LegendreBasis)
end

function SimplePolynomial(p::BasisPolynomial)
    if isa(p, SimplePolynomial)
        return p
    else
        T = eltype(p.coeffs)
        rcoeff = zero(p.coeffs)
        for i in eachindex(p.coeffs)
            pcoeff = (((basis_type(p))(i-1)).coeffs).* p.coeffs[i]
            n = length(rcoeff) - length(pcoeff)
            if n > 0
                pcoeff = vcat(pcoeff, zeros(T, n))
            end
            rcoeff += pcoeff
            
        end
        return SimplePolynomial(rcoeff)
    end
end

base_string(p::SimpleBasis) = p.order > 0 ? "x^$(p.order)" : ""
base_string(p::LegendreBasis) = "L_$(p.order)"
base_string(p::ChevyshevBase) = "T_$(p.order)"


