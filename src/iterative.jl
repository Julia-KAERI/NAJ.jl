function iteration_jacobi(
    A::AbstractMatrix, 
    b::Vector, 
    x0::Vector; 
    etol::Number = 1.0e-10,
    Maxiter::Integer = 100_000)
    n = size(A)[1]
    @assert n == size(A)[2] == size(b)[1]
    @assert Maxiter > 3
    x = zero(x0)
    
    for niter in 1:Maxiter
        @inbounds for i in 1:length(x0)
            @inbounds for j in 1:length(x0)
                if i ≠ j
                    x[i] += -A[i,j] * x0[j]
                end
            end
            @inbounds x[i] = (x[i]+b[i])/A[i, i] 
        end
        if norm(x .- x0, Inf) / norm(x, Inf) < etol
            break
        else 
            x0 = x[:]
            x = zero(x0)

        end
    end
    return x
end

function iteration_seidel(
    A::AbstractMatrix, 
    b::Vector, 
    x0::Vector; 
    etol::Number = 1.0e-10, 
    Maxiter = 100_000)
    @assert size(A)[1] == size(A)[2] == size(b)[1]
    
    x = zero(x0)
    for niter in 1:Maxiter
        @inbounds for i in 1:length(x0)
            @inbounds for j in 1:length(x0)
                if j < i
                    x[i] += - A[i, j]*x[j]
                elseif j>i
                    x[i] += - A[i, j]*x0[j] 
                end
            end
            @inbounds x[i] = (x[i]+b[i])/A[i, i] 
        end
        if norm(x .- x0, Inf) / norm(x, Inf) < etol
            break
        else 
            x0 = x[:]
            x = zero(x0)

        end
    end
    return x
end


function iteration_sor(
    A::AbstractMatrix, 
    b::Vector, 
    x0::Vector,
    ω::Real; 
    etol::Number = 1.0e-10, 
    Maxiter = 100_000)

    @assert 0 < ω < 2
    
    x = zero(x0)
    for nitter in 1:Maxiter
        @inbounds for i ∈ 1:length(x0)
            @inbounds for j ∈ 1:length(x0)
                if j < i
                    x[i] += -A[i, j] *x[j]
                elseif j > i
                    x[i] += -A[i, j] *x0[j]
                end
            end
            x[i] = (1-ω)*x0[i] + 1/A[i, i] * (ω*b[i] + ω *  x[i])
        end
        if norm(x .- x0, Inf)/norm(x, Inf)< etol    
            return x
        else 
            x0 = x[:]
            x = zero(x)
        end
    end
    return x
end


function iteration_steepest(
    A::AbstractMatrix, 
    b::Vector, 
    x0::Vector;
    etol::Number = 1.0e-5, 
    Maxiter = 100_000)

    x = similar(x0)
    for i in 1:Maxiter
        v = b - A*x0
        t = dot(v,(b-A*x0))/dot(v, (A*v))
        x = x0 + t*v
        if norm(A*x-b, Inf)<etol
            nitter = i
            println(nitter)
            return x
        else 
            x0 = x
        end
    end
    return nothing
end


function iteration_orthogonal(
    A::AbstractMatrix, 
    b::Vector, 
    x0::Vector;
    etol::Number = 1.0e-5, 
    Maxiter = 100_000)

    x = similar(x0)
    for i in 1:Maxiter
        v = b - A*x0
        t = dot(v,(b-A*x0))/dot(v, (A*v))
        x = x0 + t*v
        if norm(A*x-b, Inf)<etol
            nitter = i
            println(nitter)
            return x
        else 
            x0 = x
        end
    end
    return nothing
end