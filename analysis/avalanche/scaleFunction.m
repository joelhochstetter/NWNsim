function sf =  scaleFunction(t, T, a, gamma)
%{
    t is rescaled time 0<t<T
    Polynomial function expressing the conditions for a scaling function.
    F(0) = 0, F(1) = 0
    The basis functions are orthonormal polynomails of increasing degree
    F(t) = (a(1)*f1(t) + a(2)*f2(t) + a(3)*f3(t) + a(4)*f4(t) + a(5)*f5(t))
    
    This function returns T^gamma * F(t)
%}

    f1 = sqrt(30) .* (-1 + t).*t;
    f2 = sqrt(210) .* (-1 + t).*t.*(-1 + 2.*t);
    f3 = 3.*sqrt(10).*(-1 + t).*t.*(3 - 14.*t + 14.*t.^2);
    f4 = sqrt(2310).*(-1 + t).*t.*(-1 + 2.*t).*(1 - 6.*t + 6.*t.^2);
    f5 = 2.*sqrt(1365).*(-1 + t).*t.*(1 - 12.*t + 45.*t.^2 - 66.*t.^3 + 33.*t.^4);
    sf = (T.^gamma) .* (a(1).*f1 + a(2).*f2 + a(3).*f3 + a(4).*f4 + a(5).*f5);
    
end