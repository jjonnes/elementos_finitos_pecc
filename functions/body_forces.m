function F = body_forces(C, h, gamma, q, NoElem)
q           = quadrature(q);
Npst        = size(q, 1);
nnodes      = size(C, 1);
F           = zeros(nnodes, 1);
for i       = 1:Npst
    xi      = q (i, 1);
    eta     = q (i, 2);
    w       = q (i, 3);
    [dN, N] = quad_shape_form(NoElem, xi, eta);
    J       = C'*dN;
    F       = F + N*gamma*h*det(J)*w;
end
end
%%
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