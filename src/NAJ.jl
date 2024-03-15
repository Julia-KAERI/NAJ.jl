module NAJ

using LinearAlgebra, SparseArrays

include("polynomials/polynomials.jl")
include("interpolations.jl")
include("simples.jl")
include("calculus.jl")
include("rootfinding.jl")
include("iterative.jl")
include("ode.jl")
include("gaussian_quadrature.jl")
include("bezier.jl")

export 
    order,
    degree,
    base_type,
    base_string,
    poly_string,

    AbstractBasisPolynomial,
    SimplePolynomial,
    LegendrePolynomial,
    ChevyshevPolynomial,
    monic,
    derivative,
    integrate,
    polynomial_from_roots,
    valdermond_polynomial,
    lagrange_polynomial,
    newton_polynomial,
    least_square_poly,
    Interpolator1D,
    
    LegendrePolynomial,
    toSimplePolynomial,

    neville, 

    difference,
    integrate_trapzoidal,
    integrate_simpson_1_3,
    integrate_simpson_3_8,
    rhomberg,

    integrate_gauss_quadrature,

    rootfinding_bisection,
    rootfinding_newton,
    rootfinding_secant,
    rootfinding_regula_falci,
    
    iteration_jacobi,
    iteration_gauss_siedel,
    iteration_sor,
    iteration_steepest,
    iteration_orthogonal,

    ode_euler,
    ode_rk2,
    ode_rk4,

    Bezier

end # module NAJ