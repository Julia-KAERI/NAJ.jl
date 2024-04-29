abstract type AbstractModular{N} end
abstract type AbstractPrimeModular{N} <: AbstractModular{N} end

const max_modular = 100_000

function isprime(v::Integer)
    @assert 1<v<max_modular 
    k = floor(Int64, √v)
    result = true
    for x in 2:k
        if v % x == 0
            result = false
            break
        end
    end
    return result
end

struct PrimeMod{N} <: AbstractPrimeModular{N}
    value::Int64
    
    function PrimeMod{N}(x::T) where {N, T<:Integer}
        @assert N > 1 "Modulos must be larger than 1"
        @assert isprime(N) "Must be a prime number"
        r = (x ≥ 0) ? x%N : (x%N)+N 
        return new{N}(r)
    end
end

struct Mod{N} <: AbstractModular{N}
    val::Int64

    function Mod{N}(x::T) where {N, T<:Integer}
        @assert N > 1 "Modulos must be larger than 1"
        r = (x ≥ 0) ? x%N : (x%N)+N 
        if isprime(N)
            return PrimeMod{N}(r)
        else 
            return new{N}(r)
        end
    end
end

function Base.show(io::IO, x::AbstractModular{N}) where N 
    println(io, "$(x.value) ($N)")
end

Base.zero(a::AbstractModular{N}) where N = Mod{N}(0)
Base.one(a::AbstractModular{N}) where N = Mod{N}(1)
Base.:-(a::AbstractModular{N}) where N = Mod{N}(-a.value)

function Base.:+(a::AbstractModular{N}, b::AbstractModular{N}) where N 
    return Mod{N}(a.value+b.value)
end

function Base.:-(a::AbstractModular{N}, b::AbstractModular{N}) where N 
    return Mod{N}(a.value-b.value)
end

function Base.:*(a::AbstractPrimeModular{N}, b::AbstractPrimeModular{N}) where N 
    return Mod{N}(a.value*b.value)
end

function Base.:^(a::AbstractModular{N}, n::Integer) where N
    return Mod{N}(a.value^n)
end

"""
    minv(a::AbstractPrimeModular{N})

multiplicative inverse of a
"""
function minv(a::AbstractPrimeModular{N}) where N
    return Mod{N}(a.value^(N-2))
end

function Base.:/(a::AbstractPrimeModular{N}, b::AbstractPrimeModular{N}) where N
    if b == zero(b)
        throw(DivideError)
    end
    return a*minv(b)
end