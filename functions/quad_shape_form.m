function [dN, N] = quad_shape_form(x, xi, eta)
if x == 4
    [dN, N] = quad4_derivs (xi, eta);
elseif x == 8
    [dN, N] = quad8_derivs (xi, eta);
elseif x == 12    
    [dN, N] = quad12_derivs (xi, eta);
else
end    
end