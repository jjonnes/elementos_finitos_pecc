function K = compute_K(C, D, q, h, NoElem)
q           = quadrature(q);
Npst        = size(q, 1);
nnodes      = size(C, 1);
ndof        = 2;
K           = zeros(nnodes*ndof, nnodes*ndof);
for i       = 1:Npst
    xi      = q (i, 1);
    eta     = q (i, 2);
    w       = q (i, 3);
    B       = compute_B(C, xi, eta, NoElem);
    [dN, N] = quad_shape_form(NoElem, xi, eta);
    J       = C'*dN;
    K       = K + B'*D*B*det(J)*w*h;
end
end

function B = compute_B(C, xi, eta, NoElem)
nnodes  = size(C, 1);
ndof    = 2;
[dN, N] = quad_shape_form(NoElem, xi, eta);
J       = C'*dN;
dNdX    = dN/J;
    for i = 1: nnodes
    c = (i-1) * ndof;
    B(1, c+1) = dNdX(i,1);
    B(2, c+2) = dNdX(i,2);
    B(3, c+1) = dNdX(i,2);
    B(3, c+2) = dNdX(i,1);
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