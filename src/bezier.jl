"""
    Bezier{T}

gnenerate quadratic and cubic Bezier curve.


example
=======

```julia
c1 = Bezier([-2, 2], [2, 2], [0, 0])   # quadratic Bezier curve
c2 = Bezier([-2, 2], [2, 2], [-1, 0], [1, 0])   # cubic Bezier curve

c1(0.5)        # return bezier point for t=0.5 
```

"""
struct Bezier{T} 
    p1::Vector{T}
    p2::Vector{T}
    c1::Vector{T}
    c2::Union{Vector{T}, Nothing}

    function Bezier(
        p1::AbstractVector{T1}, 
        p2::AbstractVector{T2}, 
        c1::AbstractVector{T3}, 
        c2::Union{AbstractVector{T4}, Nothing}=nothing) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real}
        
        @assert length(p1) == length(p2) == length(c1)
        
        if c2 === nothing
            TT = promote_type(eltype(p1), eltype(p2), eltype(c1))
            return new{TT}(Vector(p1), Vector(p2), Vector(c1), nothing)
        else 
            @assert length(p1) == length(c2)
            TT = promote_type(eltype(p1), eltype(p2), eltype(c1), eltype(c2))
            return new{TT}(Vector(p1), Vector(p2), Vector(c1), Vector(c2))
        end
        
    end
end

function (b::Bezier)(t)
    if b.c2 === nothing
        return _bezier1(b.p1, b.p2, b.c1, t)
    else 
        return _bezier2(b.p1, b.p2, b.c1, b.c2, t)
    end
end

function _bezier1(p1, p2, c1, t)
    
    return (1-t)^2 .* p1 .+ (2*(1-t)*t) .*c1 .+ (t^2) .* p2
end

function _bezier2(p1, p2, c1, c2, t)
   x1 = _bezier1(p1, c2, c1, t)
   x2 = _bezier1(c1, p2, c2, t)
   return (1-t).* x1 .+ t .* x2
end
