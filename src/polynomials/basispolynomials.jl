include("coefficients.jl")

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


struct LegendreBase{T}<:PolynomialBase
    order::Integer
    coeffs::Vector{T} 

    function LegendreBase(order::Integer = 0) 
        @assert 0 ≤ order ≤ 21
        return new{Float64}(order, legendre_coefficients[order+1])
    end

    function LegendreBase{T}(order::Integer = 0) where T<:Real
        @assert 0 ≤ order ≤ 21
        return new{T}(order, legendre_coefficients[order+1])
    end
end


struct ChevyshevBase{T}<:PolynomialBase
    order::Integer
    coeffs::Vector{T}

    function ChevyshevBase(order::Integer = 0) 
        @assert 0 ≤ order ≤ 21
        return new{Float64}(order, chevyshef_coefficients[order+1])
    end

    function ChevyshevBase{T}(order::Integer = 0) where T<:Real
        @assert 0 ≤ order ≤ 21
        return new{T}(order, chevyshef_coefficients[order+1])
    end
end

struct HermiteBase{T}<:PolynomialBase
    order::Integer
    coeffs::Vector{T}

    function HermiteBase(order::Integer = 0) 
        @assert 0 ≤ order ≤ 21
        return new{Float64}(order, hermite_coefficients[order+1])
    end

    function HermiteBase{T}(order::Integer = 0) where T<:Real
        @assert 0 ≤ order ≤ 21
        return new{T}(order, hermite_coefficients[order+1])
    end
end


struct LaguerreBase{T}<:PolynomialBase
    order::Integer
    coeffs::Vector{T}

    function LaguerreBase(order::Integer = 0) 
        @assert 0 ≤ order ≤ 21
        return new{Float64}(order, laguerre_coefficients[order+1])
    end

    function LaguerreBase{T}(order::Integer = 0) where T<:Real
        @assert 0 ≤ order ≤ 21
        return new{T}(order, laguerre_coefficients[order+1])
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

@doc"""
    SimplePolynomal{T}(coeffs::AbstractVector)

Construct a polynomial from its coefficients `coeffs`, lowest order first. For example, `coeff[1]` 
be constant term and `coeff[2]` be the coefficient of `x`. `T` is the type of coefficients. The 
argument `coeffs` are convected to type `T` vector.
"""
const SimplePolynomial{T} = BasisPolynomial{T, SimpleBase}

@doc"""
    LegendrePolynomial{T}(coeffs::AbstractVector)

Construct a polynomial from Legendre polynomials. `LegendrePolynoimal([1, 2])` construct a polynmial 
function of `1 L_0(x) + 2 L_1(x)`, where `L_k(x)` is k-th order Legendre polynomial. `T` is the type 
of coefficients. The argument `coeffs` are convected to type `T` vector.
"""
const LegendrePolynomial{T} = BasisPolynomial{T, LegendreBase}

@doc"""
    ChevyshevPolynomial{T}(coeffs::AbstractVector)

Construct a polynomial from Chevyshev polynomials of type I. `ChevyshevPolynoimal([1, 2])` construct 
a polynmial function of `1 T_0(x) + 2 T_1(x)`, where `T_k(x)` is k-th order Chevyshev polynomial. 
`T` is the type of coefficients. The argument `coeffs` are convected to type `T` vector.
"""
const ChevyshevPolynomial{T} = BasisPolynomial{T, ChevyshevBase}

@doc"""
    HermitePolynomial{T}(coeffs::AbstractVector)

Construct a polynomial from Hermite polynomials (for physics). `HermitePolynoimal([1, 2])` construct
a polynmial function of `1 H_0(x) + 2 H_1(x)`, where `H_k(x)` is k-th order Hermite polynomial. 
`T` is the type of coefficients. The argument `coeffs` are convected to type `T` vector.

Note
====
Two kinds of Hermite polynomials are used(see https://en.wikipedia.org/wiki/Hermite_polynomials). This
Herimite is known for physicist's Hermite polynomials
"""
const HermitePolynomial{T} = BasisPolynomial{T, HermiteBase}

@doc"""
    LaguerrePolynomial{T}(coeffs::AbstractVector)

Construct a polynomial from Laguerre polynomials of type I. `LaguerrePolynoimal([1, 2])` construct 
a polynmial function of `1 L_0(x) + 2 L_1(x)`, where `L_k(x)` is k-th order Chevyshev polynomial. 
`T` is the type of coefficients. The argument `coeffs` are convected to type `T` vector.
"""
const LaguerrePolynomial{T} = BasisPolynomial{T, LaguerreBase}


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


function LegendreBasePolynomial(order::Integer)
    @assert order ≥ 0
    coeffs = zeros(order+1)
    coeffs[end]=1
    return LegendrePolynomial(coeffs)
end




function ChevyshevPolynomial(a::AbstractVector{T}) where T<:Real
    return BasisPolynomial(a, ChevyshevBase)
end

function ChevyshevPolynomial{P}(a::AbstractVector{T}) where {P<:Real, T<:Real}
    return BasisPolynomial{P}(a, ChevyshevBase)
end


function ChevyshevBasePolynomial(order::Integer)
    @assert order ≥ 0
    coeffs = zeros(order+1)
    coeffs[end]=1
    return ChevyshevPolynomial(coeffs)
end



function HermitePolynomial(a::AbstractVector{T}) where T<:Real
    return BasisPolynomial(a, HermiteBase)
end

function HermitePolynomial{P}(a::AbstractVector{T}) where {P<:Real, T<:Real}
    return BasisPolynomial{P}(a, HermiteBase)
end

function HermiteBasePolynomial(order::Integer)
    @assert order ≥ 0
    coeffs = zeros(order+1)
    coeffs[end]=1
    return HermitePolynomial(coeffs)
end



function LaguerrePolynomial(a::AbstractVector{T}) where T<:Real
    return BasisPolynomial(a, LaguerreBase)
end

function LaguerrePolynomial{P}(a::AbstractVector{T}) where {P<:Real, T<:Real}
    return BasisPolynomial{P}(a, LaguerreBase)
end

function LaguerreBasePolynomial(order::Integer)
    @assert order ≥ 0
    coeffs = zeros(order+1)
    coeffs[end]=1
    return LaguerrePolynomial(coeffs)
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
base_string(p::HermiteBase) = "H_$(p.order)(x)"
base_string(p::LaguerreBase) = "L_$(p.order)(x)"

poly_string(p::SimplePolynomial) = "Simple Polynomial"
poly_string(p::LegendrePolynomial) = "Legendre Polynomial"
poly_string(p::ChevyshevPolynomial) = "Chevyshev Polynomnial"
poly_string(p::HermitePolynomial) = "Hermite Polynomnial"
poly_string(p::LaguerrePolynomial) = "Laguerre Polynomnial"

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
