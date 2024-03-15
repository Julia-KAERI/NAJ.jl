
struct SimpleBase{T}<:PolynomialBase
    order::Integer
    coeffs::Vector{T} 

    function SimpleBase(order::Integer)
        @assert order ≥ 0
        r = zeros(Float64, order+1)
        r[end]=one(Float64)
        return new{Float64}(order, r)
    end

    function SimpleBase{T}(order::Integer) where T<:Real
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

struct LegendreBase{T}<:PolynomialBase
    order::Integer
    coeffs::Vector{T} 

    function LegendreBase(order::Integer = 0) 
        @assert order ≥ 0
        return new{Float64}(order, legendre_polynomial_coefficients(Float64, order))
    end

    function LegendreBase{T}(order::Integer = 0) where T<:Real
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

struct ChevyshevBase{T}<:PolynomialBase
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


const SimplePolynomial{T} = BasisPolynomial{T, SimpleBase}
const LegendrePolynomial{T} = BasisPolynomial{T, LegendreBase}
const ChevyshevPolynomial{T} = BasisPolynomial{T, ChevyshevBase}


function SimplePolynomial(a::AbstractVector{T}) where {T<:Real}
    return BasisPolynomial(a, SimpleBase)
end

function SimplePolynomial{P}(a::AbstractVector{T}) where {P<:Real, T<:Real}
    return BasisPolynoimal{P}(a, SimpleBase)
end

function LegendrePolynomial(a::AbstractVector{T}) where T<:Real
    return BasisPolynomial(a, LegendreBase)
end

function LegendrePolynomial{P}(a::AbstractVector{T}) where {P<:Real, T<:Real}
    return BasisPolynomial{P}(a, LegendreBase)
end

function ChevyshevPolynomial(a::AbstractVector{T}) where T<:Real
    return BasisPolynomial(a, ChevyshevBase)
end


function ChevyshevPolynomial{P}(a::AbstractVector{T}) where {P<:Real, T<:Real}
    return BasisPolynomial{P}(a, ChevyshevBase)
end

function SimplePolynomial(p::BasisPolynomial)
    if isa(p, SimplePolynomial)
        return p
    else
        T = eltype(p.coeffs)
        rcoeff = zero(p.coeffs)
        for i in eachindex(p.coeffs)
            pcoeff = (((base_type(p))(i-1)).coeffs).* p.coeffs[i]
            n = length(rcoeff) - length(pcoeff)
            if n > 0
                pcoeff = vcat(pcoeff, zeros(T, n))
            end
            rcoeff += pcoeff
            
        end
        return SimplePolynomial(rcoeff)
    end
end

base_string(p::SimpleBase) = p.order > 0 ? "x^$(p.order)" : ""
base_string(p::LegendreBase) = "P_$(p.order)(x)"
base_string(p::ChevyshevBase) = "T_$(p.order)(x)"

poly_string(p::SimplePolynomial) = "Simple Polynomial"
poly_string(p::LegendrePolynomial) = "Legendre Polynomial"
poly_string(p::ChevyshevPolynomial) = "Chevyshev Polynomnial"

function Base.show(io::IO, p::AbstractBasisPolynomial{T, X}) where {T, X}
    result = ""
    n = length(p)
    if n == 1 && iszero(p.coeffs[1])
        result = "0"
    else 
        isfirst = true
        for (i, v) in enumerate(p.coeffs[1:end])
            vp = string(abs(v))
            if v > zero(T) && ~(isfirst)
                result *= " + " * vp * " " * base_string(base_type(p)(i-1))
            elseif v< zero(T) && ~(isfirst)
                result *= " - " * vp * " " * base_string(base_type(p)(i-1))
            elseif v > zero(T)
                result = vp * " " * base_string(base_type(p)(i-1))
                isfirst = false
            elseif v<zero(T)
                result = " -"*vp * " " * base_string(base_type(p)(i-1))
                isfirst = false
            end

        end
    end
    println(io, poly_string(p)*"{$(eltype(p.coeffs))}($(result[1:end]))")
end
