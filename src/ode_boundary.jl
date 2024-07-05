


function linear_shooting(p::Function, q::Function, r::Function, boundary, condition, Npoints = 100)
    
    @assert length(boundary) == length(condition) == 2
    
    de1(t, x) = [0 1 0 ; q(t) p(t) r(t);  0 0 0] * x
    de2(t, x) = [0 1; q(t) p(t)] * x

    t = range(boundary[1], boundary[2], length = Npoints)

    α, β = condition[1:2]

    tn, x1 = ode_rk4(de1, boundary[1], [α, 0, 1], Npoints-1, step(t))
    tn, x2 = ode_rk4(de2, boundary[1], [0, 1], Npoints-1, step(t))
    
    r = @. x1[1, :] + (β-x1[1, end])/(x2[1, end]) * x2[1, :]

    return tn, r
end

