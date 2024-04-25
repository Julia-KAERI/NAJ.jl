struct LogSpacingRange{T}
    start::T
    stop::T
    base::T
    size::Int64

    function LogSpacingRange(start::Real, stop::Real, size::Int64, base::Real)
        @assert base > 1
        @assert size ≥ 1
        ftype = promote_type(typeof(start), typeof(stop))
        if ~(ftype <: AbstractFloat)
            ftype = Float64
        end

        return new{ftype}(start, stop, base, size)
    end
end

Base.size(p::LogSpacingRange) = (p.size, )
Base.length(p::LogSpacingRange) = p.size
Base.iterate(p::LogSpacingRange, state=1) = state > p.size ? nothing : (p[state], state+1)
Base.eltype(p::LogSpacingRange{T}) where {T} = T
Base.IteratorSize(::LogSpacingRange{T}) where {T} = Base.HasLength()
Base.IteratorEltype(::LogSpacingRange{T}) where {T} = Base.HasEltype()
Base.firstindex(p::LogSpacingRange) = (p.base)^(p.start)
Base.lastindex(p::LogSpacingRange) = (p.base)^(p.stop)

function Base.getindex(p::LogSpacingRange, i::Int64) 
    @assert 0 < i ≤ p.size
    if p.size == 1
        return p.start
    else 
        r = p.start + (p.stop-p.start)/(p.size-1)*(i-1)
        return (p.base)^r
    end
end

function logspace(a, b, n::Integer=10, base::Real=10)
    return LogSpacingRange(a, b, n, base)
end


"""
        logspace(a, b, n::Integer=10, base::Real=10)
    
Generate ``n`` element iterator for ``base^a`` to ``base^b``.

# Examples

```jldoctest
In [4]: for v in logspace(-1, 3, 5)
        @show v
        end
v = 0.1
v = 1.0
v = 10.0
v = 100.0
v = 1000.0
```
"""
function logspace(a, b, n::Integer=10, base::Real=10)
    return LogSpacingRange(a, b, n, base)
end