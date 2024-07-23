function [dN, N] = lin4_derivs (r)
syms xi eta
n = [1.0/16.0*( -9.0*xi^3 + 9.0*xi*xi +      xi - 1.0)
     1.0/16.0*(  9.0*xi^3 + 9.0*xi*xi -      xi - 1.0)
     1.0/16.0*( 27.0*xi^3 - 9.0*xi*xi - 27.0*xi + 9.0)
     1.0/16.0*(-27.0*xi^3 - 9.0*xi*xi + 27.0*xi + 9.0)];
dN = subs([diff(n, xi)], xi, r);
N  = subs(n, xi, r);
end