function difference(f::Function, x::T, ϵ = 1.0e-6, npts = 3) where T<:Real
    @assert npts ∈ (3, 5)
    eps = convert(T, ϵ)
    if npts == 3
        r = (f(x+eps) - f(x-eps))/(2*eps)
    elseif npts == 5
        r = (f(x-2*eps) - 8*f(x-eps) + 8*f(x+eps)-f(x+2*eps))/(12*eps)
    else 
        error("npts be 3 or 5")
    end
    return r
end

function integrate_trapzoidal(
    f::Function, 
    a::Number, 
    b::Number, 
    n::Integer)
    
    a, b = minmax(a, b)
    h = (b-a)/(n-1)
    x = range(a, b, length = n)
    ff = f.(x)
    result = 0.5*(ff[1]+ff[end])
    result += sum(ff[2:end-1])
    return result*h
end

function integrate_simpson_1_3(
    f::Function, 
    a::Number, 
    b::Number, 
    n::Integer)
    
    @assert n %2 == 1
    a, b = minmax(a, b)
    h = (b-a)/(n-1)
    x = range(a, b, length = n)
    ff = f.(x)
    result = (ff[1]+ff[end])
    result += 4*sum(ff[2:2:end-1])
    result += 2*sum(ff[3:2:end-2])
    return result * h/3
end

function integrate_simpson_3_8(
    f::Function, 
    a::Number, 
    b::Number, 
    n::Integer)
    
    @assert n %3 == 1
    a, b = minmax(a, b)
    h = (b-a)/(n-1)
    x = range(a, b, length = n)
    println("h=$h, dx=$(x[10]-x[9])")
    ff = f.(x)
    result = (ff[1]+ff[end])
    result += 3*sum(ff[2:3:end-2])
    result += 3*sum(ff[3:3:end-1])
    result += 2*sum(ff[4:3:end-3])
    return result * h * 3 / 8
end

function rhomberg(f::Function, a::Number, b::Number, order::Integer = 4)
    @assert order > 2
    N = zeros(order, order) 
    N[:, 1]= [integrate_trapzoidal(f, a, b, 2^k) for k in 1:order]
    for m = 1:order, k = 1:order-m
        N[k, m+1] = 1.0/(4^m-1.0) *(4^m * N[k+1, m]- N[k, m]) 
    end
    return N[1, end]
end