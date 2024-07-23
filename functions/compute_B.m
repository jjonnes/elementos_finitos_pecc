function B = compute_B(C, q, NoElem, ndof)
q           = quadrature(q);
Npst        = size(q, 1);
nnodes      = size(C, 1);
for i       = 1:Npst
    xi      = q (i, 1);
    eta     = q (i, 2);
    w       = q (i, 3);
    [dN, N] = quad_shape_form(NoElem, xi, eta);
    J       = C'*dN;
    dNdX    = dN/J;
    for j = 1: nnodes
    c = (j-1) * ndof;
    B(1, c+1) = dNdX(j,1);
    B(2, c+2) = dNdX(j,2);
    B(3, c+1) = dNdX(j,2);
    B(3, c+2) = dNdX(j,1);
    end
end
end

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