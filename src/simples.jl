function neville(
    xp::AbstractVector{T1}, 
    yp::AbstractVector{T2}, 
    x::Real) where {T1<:Real, T2<:Real}
    @assert length(xp) == length(yp)
    P = copy(yp)
    for i in 1:(length(xp)-1)
        P = [((x-xp[j+i])*P[j] - (x-xp[j])*P[j+1])/(xp[j]-xp[j+i]) for j in 1:(length(xp)-i)]
    end
    return P[1]
end