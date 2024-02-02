module NAJ

include("simplepolynomials.jl")
include("interpolations.jl")

export 
    SimplePolynomial,
    monic,
    derivative,
    integrate,
    polynomial_from_roots,
    valdermond_polynomial,
    lagrange_polynomial,
    newton_polynomial,
    Interpolator1D

end # module NAJ