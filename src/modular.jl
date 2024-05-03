abstract type AbstractMod{T, N} end 

const max_modular = 1_000_000_000

function isprime(v::Integer)
    @assert 1<v<max_modular 
    result = true
    for x in 2:floor(Int64, √v)
        if v % x == 0
            result = false
            break
        end
    end
    return result
end

struct PrimeMod{T, N} <: AbstractMod{T, N}
    value::T
    
    function PrimeMod{T, N}(x::Integer) where {T<:Integer, N}
        @assert N > 1 "Modulos must be larger than 1"
        @assert isprime(N) "Must be a prime number"
        r = (x ≥ 0) ? x%N : (x%N)+N 
        return new{T, N}(r)
    end
end

struct Mod{T, N} <: AbstractMod{T, N}
    value::T

    function Mod{T, N}(x::Integer) where {T<:Integer, N}
        @assert N > 1 "Modulos must be larger than 1"
        x = T(x)
        r = (x ≥ 0) ? x%N : (x%N)+N 
        if isprime(N)
            return PrimeMod{T, N}(r)
        else 
            return new{T, N}(r)
        end
    end
end

Mod{N}(x::Integer) where {N} = Mod{Int64, N}(x)
PrimeMod{N}(x::Integer) where {N} = PrimeMod{Int64, N}(x)


function Base.show(io::IO, x::AbstractMod{T, N}) where {T, N} 
    println(io, "$(x.value)_$N")
end

Base.zero(a::AbstractMod{T, N}) where {T, N}  = Mod{T, N}(0)
Base.one(a::AbstractMod{T, N}) where {T, N}  = Mod{T, N}(1)
Base.:-(a::AbstractMod{T, N}) where {T, N}  = Mod{T, N}(-a.value)
Base.:isequal(a::AbstractMod{T, N}, b::AbstractMod{T, N}) where {T, N} = (a.value == b.value)

function Base.:+(a::AbstractMod{T, N}, b::AbstractMod{T, N}) where {T, N} 
    return Mod{T, N}(a.value+b.value)
end

function Base.:-(a::AbstractMod{T, N}, b::AbstractMod{T, N}) where {T, N}  
    return Mod{T, N}(a.value-b.value)
end

function Base.:*(a::AbstractMod{T, N}, b::AbstractMod{T, N}) where {T, N}  
    return Mod{T, N}(a.value*b.value)
end

function Base.:inv(a::PrimeMod{T, N}) where {T, N} 
    return Mod{T, N}(a.value^(N-2))
end

function Base.:^(a::AbstractMod{T, N}, n::Integer) where {T, N} 
     return Mod{T, N}(a.value^n)
end

function Base.:/(a::PrimeMod{T, N}, b::PrimeMod{T, N}) where {T, N} 
    if b == zero(b)
        throw(DivideError)
    end
    return a*inv(b)
end