%Exercises 5.2
%1.The strain in a 1D truss element is given by ε = ∂u/∂x. Mountthe corresponding matrix B in terms of nodal coordinates x1 and x2.
clear all
clc

syms x1 x2 xi

C = [x1; x2];

B = simplify(compute_B(C, xi));

function B = compute_B(C, xi)
nnodes = size(C, 1);
ndof = 1;
dN   = lin2_derivs(xi);
J    = C'*dN;
dNdX = dN*inv(J);
    for i = 1: nnodes
    c = (i-1) * ndof;
    B(1, c+1) = dNdX(i,1);
    end
end

function dN = lin2_derivs (xi)
n = [ 1/2 - xi/2
      1/2 + xi/2];
dN = [diff(n, xi)];
end
%Exercises 5.2
%2. Mount the matrix B for a three-node bar element given thenodal coordinates x1, x2 and x3.
clear all
clc

syms x1 x2 x3 xi

C = [x1; x2; x3];

B = compute_B(C, xi);

function B = compute_B(C, xi)
nnodes = size(C, 1);
ndof = 1;
dN   = lin3_derivs(xi);
J    = C'*dN;
dNdX = dN*inv(J);
    for i = 1: nnodes
    c = (i-1) * ndof;
    B(1, c+1) = dNdX(i,1);
    end
end

function dN = lin3_derivs (xi)
n = [ xi*xi/2 - xi/2
      xi*xi/2 + xi/2
            1 - xi*xi];
dN = [diff(n, xi)];
end
%Exercises 5.2
%3. Find an expression for the matrix B of a three-node bar element located in 2D space as a function of the Jacobian and the local coordinate ξ. Consider ε = ∂u/∂s wheres is the curviline arco ordinate along the bar path.
clear all
clc

syms x1 x2 x3 
syms y1 y2 y3
syms xi

C = [x1 y1; x2 y2; x3 y3];

B = compute_B(C, xi);

function B = compute_B(C, xi)
nnodes = size(C, 1);
ndof = 2;
dN   = lin3_derivs(xi);
J    = (C'*dN);
dNdX = dN * (J/(norm(J)*norm(J)))';
    for i = 1: nnodes
    c = (i-1) * ndof;
    B(1, c+1) = dNdX(i,1);
    B(2, c+2) = dNdX(i,2);
    end
end

function dN = lin3_derivs (xi)
n = [ xi*xi/2 - xi/2
      xi*xi/2 + xi/2
            1 - xi*xi];
dN = [diff(n, xi)];
end
%Exercises 5.2
%4.Find the matrix B for a three-node triangular element given nodal coordinates(x1,y1),(x2,y2) and (x3, y3). Why matrix B is constant along the element and what does it imply?
clear all
clc

syms x1 x2 x3 
syms y1 y2 y3
syms xi eta

C = [x1 y1; x2 y2; x3 y3];

B = compute_B(C, xi, eta);

function B = compute_B(C, xi, eta)
nnodes = size(C, 1);
ndof = 2;
dN   = tri3_derivs(xi, eta);
J    = C'*dN;
dNdX = dN/J;
    for i = 1: nnodes
    c = (i-1) * ndof;
    B(1, c+1) = dNdX(i,1);
    B(2, c+2) = dNdX(i,2);
    B(3, c+1) = dNdX(i,2);
    B(3, c+2) = dNdX(i,1);
    end
end

function dN = tri3_derivs (xi, eta)
n = [ 1 - xi - eta
          xi
               eta];
dN = [diff(n, xi) diff(n, eta)];
end
%Exercises 5.2
%5. For the quadrilateral element shown below, given the vector of nodal displacements U, find the strain components at point(ξ,η) = (−1√3,1√3). U = [0.0    0.0    0.01    0.01    0.015    0.015    0.0    0.015]
clear all
clc

syms xi eta epsilon
x = -1/sqrt(3);
y = 1/sqrt(3);
C = [0 0; 4 2; 4 4; 0 2];
U = [0.0  0.0  0.01  0.01  0.015  0.015  0.0  0.015];
b = compute_B(C, xi, eta);

B = subs(b, [xi eta], [x y]);
d =  B * U';
D = vpa(d, 5);

function B = compute_B(C, xi, eta)
nnodes = size(C, 1);
ndof = 2;
dN   = quad4_derivs(xi, eta);
J    = C'*dN;
dNdX = dN/J;
    for i = 1: nnodes
    c = (i-1) * ndof;
    B(1, c+1) = dNdX(i,1);
    B(2, c+2) = dNdX(i,2);
    B(3, c+1) = dNdX(i,2);
    B(3, c+2) = dNdX(i,1);
    end
end

function dN = quad4_derivs (xi, eta)
n = [1.0/4.0 * (1 - xi) * (1 - eta)
     1.0/4.0 * (1 + xi) * (1 - eta)
     1.0/4.0 * (1 + xi) * (1 + eta)
     1.0/4.0 * (1 - xi) * (1 + eta)];
dN = [diff(n, xi) diff(n, eta)];
end

