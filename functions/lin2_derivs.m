function [dN, N] = lin2_derivs (r)

syms xi eta
n = [1.0/2.0*(1 - xi)
     1.0/2.0*(1 + xi)];
dN = subs([diff(n, xi)], xi, r);
N  = subs(n, xi, r);
end