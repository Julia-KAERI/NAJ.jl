# include("simplepolynomials.jl")


function get_cubic_spline_coefficients(xp, yp, bc_kind = :cubic_spline_naturalbc, bc=nothing)
    @assert length(xp) == length(yp)
    @assert bc_kind ∈ (:cubic_spline_clampedbc, :cubic_spline_naturalbc)
    @assert (bc === nothing) || (length(bc) == 2)
    T = eltype(xp)
    N = length(xp)
    M = spzeros(T, (4*N-4, 4*N-4))
    Y = zeros(T, (4*N-4, 1))
    for i in 1:N-2 
        M[4*(i-1)+1, (4*(i-1)+1):4*i] = [one(T) xp[i] (xp[i])^2 (xp[i])^3]
        M[4*(i-1)+2, 4*(i-1)+1:4*i] = [one(T) xp[i+1] (xp[i+1])^2 (xp[i+1])^3]
        M[4*(i-1)+3, 4*(i-1)+2:4*(i+1)] = [1 2*xp[i+1] 3*(xp[i+1])^2 0 -1 -2*xp[i+1] -3*(xp[i+1])^2]
        M[4*(i-1)+4, 4*(i-1)+3:4*(i+1)] = [2 6*(xp[i+1]) 0 0 -2 -6*(xp[i+1])]
        
        Y[4*(i-1)+1] = yp[i]
        Y[4*(i-1)+2] = yp[i+1]
        
        if i == N-2
            println(4*N-4, ", ", 4*(i-1)+4)
        end
    end
    
    M[end-3, end-3:end] = [one(T) xp[end-1] (xp[end-1])^2 (xp[end-1])^3] 
    M[end-2, end-3:end] = [one(T) xp[end] (xp[end])^2 (xp[end])^3] 
    
    if bc_kind == :cubic_spline_naturalbc
        
        M[end-1, 3:4] = [2 6*xp[1]]
        M[end, end-1:end] = [2 6*xp[end]]
        Y[end-3] = yp[end-1]
        Y[end-2]=  yp[end]
    else
        M[end-3, end-3:end] = [one(T) xp[end-1] (xp[end-1])^2 (xp[end-1])^3] 
        M[end-2, end-3:end] = [one(T) xp[end] (xp[end])^2 (xp[end])^3] 
        M[end-1, 2:4] = [1 2*xp[1] 3(xp[1])^2]
        M[end, end-2:end] = [1 2*xp[end] 3(xp[end])^2]
        Y[end-3] = yp[end-1]
        Y[end-2]=  yp[end]
        Y[end-1] = bc[1]
        Y[end] = bc[2]
    end

    return (M\Y)[:,1]

end

mutable struct Interpolator1D{T}
    xp :: Vector{T}
    yp :: Vector{T}
    kind :: Symbol
    bc :: Union{Vector{T}, Nothing}
    coeffs ::Union{Vector{T}, Nothing}

    function Interpolator1D(
        xp::AbstractVector{T}, 
        yp::Vector{S}, 
        kind::Symbol, 
        bc::Union{Nothing, Vector}=nothing, 
        ) where {T <: Real, S<: Real}
        @assert kind ∈ (:nearest, :linear, :cubic, :cubic_spline_naturalbc, :cubic_spline_clampedbc)
        @assert length(xp) == length(yp)
        
        N = promote_type(T, S)
        xp = convert.(N, xp)
        yp = convert.(N, yp)

        if kind ∈ (:cubic_spline_naturalbc, :cubic_spline_clampedbc)
            coeffs = get_cubic_spline_coefficients(xp, yp, kind, bc)
        else 
            bc = nothing
            coeffs = nothing
        end
        new{N}(Vector(xp), Vector(yp), kind, bc, coeffs)
    end
end

function (p::Interpolator1D)(x::T) where T <: Number
    if x < p.xp[1] || x > p.xp[end]
        result = zero(p.yp[1])
    else
        ind = findfirst(xs->(xs>=x), p.xp)-1
        
        if p.kind == :nearest
            if ind == 0
                result = p.yp[1]
            elseif x-p.xp[ind] > p.xp[ind+1]-x
                result = p.yp[ind+1]
            else 
                result = p.yp[ind]
            end
        elseif p.kind == :linear
            if ind == 0
                result = p.yp[1]
            else
                result = (p.yp[ind+1]-p.yp[ind])/(p.xp[ind+1] - p.xp[ind])*(x - p.xp[ind]) + p.yp[ind]
            end
        elseif p.kind == :cubic
            N = length(xp) 
            if ind ∈ (0, 1)
                xs, ys = p.xp[1:4], p.yp[1:4]
            elseif ind ∈ (N, N-1)
                xs, ys = p.xp[end-3:end], p.yp[end-3:end]
            else 
                xs, ys = p.xp[ind-1:ind+2], p.yp[ind-1:ind+2]
            end
            result = newton_polynomial(xs, ys)(x)

        elseif p.kind ∈ (:cubic_spline_naturalbc, :cubic_spline_clampedbc)
            if ind == 0
                result = yp[1]
            else 
                result = evalpoly(x, p.coeffs[4*ind-3:4*ind])
            end
        else 
            result = zero(T)
        end            
    end
    return result
end
