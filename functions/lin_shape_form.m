function [dN, N] = lin_shape_form(x, xi)
if x == 2
    [dN, N] = lin2_derivs (xi);
elseif x == 3
    [dN, N] = lin3_derivs (xi);
elseif x == 4    
    [dN, N] = lin4_derivs (xi);
else
end    
end