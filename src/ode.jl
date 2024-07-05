
function ode_euler(
    fp, 
    t1::Real, 
    x1::Vector{<:Real}, 
    Npoints::Integer, 
    h = 1.0e-6)
    
    tn = t1 .+ collect(0:1:(Npoints-1)) * h
    xn = zeros((length(x1), length(tn)))
    xn[:,1] = x1
    for i in 1:(Npoints-1)
        @inbounds xn[:, i+1] = xn[:, i] .+ fp(tn[i], xn[:, i]).*h
    end
    return tn, xn
end

# function ode_euler(
#     fp,
#     t1::Real, 
#     x1::Real,
#     Npoints::Integer, 
#     h = 1.0e-6)

#     tn, xn =  ode_euler(fp, t1, [x1], Npoints, h)
#     return tn, xn[1,:]
# end

function ode_rk2(
    fp, 
    t1::Real, 
    x1::Vector{<:Real}, 
    Npoints::Integer, 
    h = 1.0e-6) 
    tn = t1 .+ collect(0:1:(Npoints-1)) * h
    xn = zeros((length(x1), length(tn)))
    xn[:,1] = x1
    for i in 1:(Npoints-1)
        @inbounds k1 = fp(tn[i], xn[:,i])
        @inbounds k2 = fp(tn[i] + h/2, xn[:, i] .+ k1.*(h/2))
        @inbounds xn[:, i+1] = xn[:,i] .+ (k1 .+ k2) .*(h/2)
    end 
    return tn, xn
end

function ode_rk4(
    fp,
    t1::Real, 
    x1::Vector{<:Real}, 
    Npoints::Integer, 
    h = 1.0e-6) 
    tn = t1 .+ collect(0:1:(Npoints-1)) * h
    xn = zeros((length(x1), length(tn)))
    xn[:,1] = x1
    for i in 1:(Npoints-1)
        @inbounds k1 = fp(tn[i], xn[:, i])
        @inbounds k2 = fp(tn[i] + h/2, xn[:, i] .+ k1.*(h/2))
        @inbounds k3 = fp(tn[i] + h/2, xn[:, i] .+ k2 .*(h/2))
        @inbounds k4 = fp(tn[i] + h, xn[:, i] .+ k3 .* h)
        @inbounds xn[:, i+1] = xn[:, i] .+ (k1 .+ (2.0 .* k2) .+ (2.0 .* k3) .+ k4).*(h/6)
    end
    return tn, xn
end

ode_initial_methods = Dict(:euler => ode_euler, :rk2 => ode_rk2, :rk4 => ode_rk4)

mutable struct InitialValueOdeProblem
    diffeq::Function
    t1::Number
    Npoints::Integer
    h::Real
    initial_value::Vector

    function InitialValueOdeProblem(diffeq::Function, t1::Number, initial_value,  Npoints::Integer, h::Real)
        @assert Npoints > 3
        @assert h > 0
        if isa(initial_value, Real) 
            initial_value = [initial_value, ]
        end
        return new(diffeq, t1, Npoints, h, initial_value)
    end
end

function solve(ode::InitialValueOdeProblem, method::Symbol)
    @assert method ∈ keys(ode_initial_methods)
    return ode_initial_methods[method](ode.diffeq, ode.t1, ode.initial_value, ode.Npoints, ode.h)

end

