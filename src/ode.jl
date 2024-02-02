function ode_euler(
    fp::Function, 
    t1::Real, 
    x1::Vector{<:Real}, 
    Npoints::Integer, 
    h = 1.0e-6)
    
    tn = t1 .+ collect(0:1:(Npoints-1)) * h
    xn = zeros((length(x1), length(tn)))
    xn[:,1] = x1
    for i in 1:(Npoints-1)
        xn[:, i+1] = xn[:, i] .+ fp(tn[i], xn[:, i]).*h
    end
    return tn, xn
end

function ode_euler(
    fp::Function,
    t1::Real, 
    x1::Real,
    Npoints::Integer, 
    h = 1.0e-6)

    tn, xn =  ode_euler(fp, t1, [x1], Npoints, h)
    return tn, xn[1,:]
end

function ode_rk2(f::Function, 
    t1::Real, 
    x1::Vector{<:Real}, 
    Npoints::Integer, 
    h = 1.0e-6) 
    tn = t1 .+ collect(0:1:(Npoints-1)) * h
    xn = zeros((length(x1), length(tn)))
    xn[:,1] = x1
    for i in 1:(Npoints-1)
        k1 = f(tn[i], xn[:,i])
        k2 = f(tn[i] + h/2, xn[:, i] .+ k1.*(h/2))
        xn[:, i+1] = xn[:,i] .+ (k1 .+ k2) .*(h/2)
    end 
    return tn, xn
end

function ode_rk4(f::Function, 
    t1::Real, 
    x1::Vector{<:Real}, 
    Npoints::Integer, 
    h = 1.0e-6) 
    tn = t1 .+ collect(0:1:(Npoints-1)) * h
    xn = zeros((length(x1), length(tn)))
    xn[:,1] = x1
    for i in 1:(Npoints-1)
        k1 = f(tn[i], xn[:, i])
        k2 = f(tn[i] + h/2, xn[:, i] .+ k1.*(h/2))
        k3 = f(tn[i] + h/2, xn[:,i] .+ k2 .*(h/2))
        k4 = f(tn[i] + h, xn[:, i] .+ k3 .* h)
        xn[:, i+1] = xn[:, i] .+ (k1 .+ (2.0 .* k2) .+ (2.0 .* k3) .+ k4).*(h/6)
    end
    return tn, xn
end