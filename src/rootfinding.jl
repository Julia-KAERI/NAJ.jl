function rootfinding_bisection(
    f,                              # 함수
    a::Real,                        # 구간값 1
    b::Real,                        # 구간값 2
    xtol::Real = 1.0e-8,            # 해의 오차의 허용 범위
    etol::Real = 1.0e-10)           # 해의 함수값의 허용 범위
    Niter = 0
    a, b = minmax(a, b)
    f(a)*f(b) <= 0 || error("f(a)*f(b) should be negative") 
    c = (a+b)/2
    while ((b-a) > 2*xtol) || (abs(f(c))<etol)
        Niter +=1
        
        if f(c) == 0.0
            break
        elseif f(a)*f(c) < 0 
            a, b = a, c
        else 
            a, b = c, b 
        end
        c = (a+b)/2
    end
    return c, Niter
end

function rootfinding_newton(
    f,                              # 함수
    df::Function,                   # 도함수
    p::Real,                        # 시작값
    MaxIter::Int64=100_000,         # 최대 반복 횟수
    etol::Real = 1.0e-8,            # 해의 함수값의 허용 범위
    dfmin::Real = 1.0e-6)           # 미분값의 절대값의 허용되는 최소값
    
    Niter = 0
    for i in 1:MaxIter
        if abs(f(p)) < etol
            break
        elseif abs(df(p)) < dfmin 
            error("df ≈ 0.0")
        end
        p = p - f(p)/df(p)
        Niter += 1       
        
        if abs(f(p)) < etol 
            return (p, Niter)
        end
    end
    error("최대 반복 횟수 $MaxIter 에 도달하였으나 답을 찾지 못함")
end

function rootfinding_secant(
    f,                          # 함수
    p0::Real,                   # 시작값 1
    p1::Real,                   # 시작값 2
    MaxIter::Int64 = 100_000,   # 최대 반복 횟수
    etol::Real = 1.0e-8,        # 해의 함수값의 허용 범위
    dfmin::Real = 1.0e-6)       # 도함수의 근사값의 절대값에 허용되는 최소값
    
    Niter = 0
    
    for i in 1:MaxIter
        gx = (f(p1)-f(p0))/(p1-p0)
        if abs(f(p1)) < etol 
            return p1, Niter
        elseif abs(gx)<dfmin
            error("df ≈ 0.0")
        end
        p0, p1 = p1, p1 - f(p1)/gx
        Niter += 1
    end
    error("최대 반복 횟수 $MaxIter 에 도달하였으나 답을 찾지 못함.")
end

function rootfinding_regula_falci(
    f,                              # 함수
    a::Real,                        # 구간값 1
    b::Real,                        # 구간값 2
    MaxIter::Int64 = 100_000,       # 최대 반복 회수
    xtol::Real = 1.0e-8,            # 해의 오차의 허용 범위
    etol::Real = 1.0e-8,            # 해의 함수값의 허용 범위
    dfmin::Real = 1.0e-6)           # 도함수의 근사값의 절대값에 허용되는 최소값

    a, b = minmax(a, b)
    @assert f(a)*f(b) < 0

    Niter = 0

    for i in 1:MaxIter
        Niter +=1
        gx =  (f(b)-f(a))/(b-a)
        c = b-f(b)/gx
        if (abs(b-a)<xtol) || (abs(f(c))<etol)
            return c, Niter
        elseif abs(gx) < dfmin 
            error("df ≈ 0.0")
        end
        
        
        if f(b)*f(c) < 0 
            a, b = b, c
        else 
            a, b = a, c
        end
    end
    error("최대 반복 횟수 $MaxIter 에 도달하였으나 답을 찾지 못함.")
end