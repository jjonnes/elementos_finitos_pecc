function [dN, N] = lin3_derivs (r)

syms xi eta
n = [1.0/2.0*(xi*xi - xi)
     1.0/2.0*(xi*xi + xi)
     1.0 -         xi*xi];
dN = subs([diff(n, xi)], xi, r);
N  = subs(n, xi, r);
end