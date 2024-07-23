%% Exercises 4.4
%1. For a three-node bar element with coordinates(x1,y1) = (1,1),(x2,y2) =(3,1)and(x3,y3) = (2,2), determine an expression for its orientationas a function of the local coordinateξ. Find the orientation vector at node 3.
clear all
clc

n = 3;
xi = linspace (-1, 1, n);

C = [1 1; 3 1; 2 2];
for i = 1:n
    J(:, :, i) = C'* lin_deriv (xi(i));
    Jacobian(:, :, i) = norm (J(:, :, i));
    r(:, :, i) =  J(:, :, i)/Jacobian(:, :, i);
end
r = r(:, :, 2)';

function dn = lin_deriv (xi)
dn = [xi - 1/2
      xi + 1/2
      -2*xi];
end
%% Exercises 4.4
%2. Using integration, calculate the length of the element in the last exer-cise.
clear all
clc

syms xi 

C = [1 1; 3 1; 2 2];
J = C'* lin_deriv (xi);
eqs = norm(J);
L = int(eqs, -1, 1);

function dn = lin_deriv (xi)
dn = [xi - 1/2
      xi + 1/2
      -2*xi];
end
%% Exercises 4.4
%3. Study the limitations for the location of intermediate nodes in a four-node bar element. Assume that the nodes are placed simmetrically.
clear all
clc

syms xi alpha L

C = [-L/2 ; L/2 ; -L/6; alpha*L];
j = C'* lin_deriv (xi);
simplify(j);
J = subs(j, L, 1) == 0;
alfa = [subs(J, xi, -1/3); ...
        subs(J, xi,  1  )];
for i = 1:2
    A(i, 1) = solve (alfa(i));
end
fprintf ('[%0.2f , %0.2f]*L\n',A)

function dn = lin_deriv (xi)
n = [1.0/16.0*( -9.0*xi^3 + 9.0*xi*xi +      xi - 1.0)
     1.0/16.0*(  9.0*xi^3 + 9.0*xi*xi -      xi - 1.0)
     1.0/16.0*( 27.0*xi^3 - 9.0*xi*xi - 27.0*xi + 9.0)
     1.0/16.0*(-27.0*xi^3 - 9.0*xi*xi + 27.0*xi + 9.0)];
dn = diff(n);
end
%% Exercises 4.4
%4. By inspection determine the Jacobian norm of the four-node elementshown below. Then, compare the result with the one obtained by usingthe Jacobian matrix.

clear all
clc

syms xi eta
n = 2;
x = linspace (-1, 1, n);

C = [0 0 0; 1 0 0; 1 1 1; 0 1 1];
dn = quad4_deriv (xi, eta);
J1 = C' * subs (dn, [xi eta], [x(1) x(1)]);
J2 = C' * subs (dn, [xi eta], [x(2) x(1)]);
J3 = C' * subs (dn, [xi eta], [x(2) x(2)]);
J4 = C' * subs (dn, [xi eta], [x(1) x(2)]);

NormJ1 = sqrt((det(J1([1 2], [1 2])))^2 + ...
              (det(J1([2 3], [1 2])))^2 + ...
              (det(J1([3 1], [1 2])))^2);
NormJ2 = sqrt((det(J2([1 2], [1 2]))^2) + ...
              (det(J2([2 3], [1 2]))^2) + ...
              (det(J2([3 1], [1 2]))^2));
NormJ3 = sqrt((det(J3([1 2], [1 2]))^2) + ...
              (det(J3([2 3], [1 2]))^2) + ...
              (det(J3([3 1], [1 2]))^2));
NormJ4 = sqrt((det(J4([1 2], [1 2]))^2) + ...
              (det(J4([2 3], [1 2]))^2) + ...
              (det(J4([3 1], [1 2]))^2));

function dn = quad4_deriv (xi, eta)
n = [1.0/4.0 * (1 - xi) * (1 - eta)
     1.0/4.0 * (1 + xi) * (1 - eta)
     1.0/4.0 * (1 + xi) * (1 + eta)
     1.0/4.0 * (1 - xi) * (1 + eta)];
dn = [diff(n, xi), diff(n, eta)];
end
%% Exercises 4.4
%5. Compute the Jacobian norm as function of ξ and η of the element in the last exercise changing the coordinates of node 3 to x3= (1,1,0).Later calculate the surface area by integration.

clear all
clc

syms xi eta

C = [0 0 0; 1 0 0; 1 1 0; 0 1 1];
dn = quad4_deriv (xi, eta);
J = C' * dn;

NormJ = sqrt((det(J([1 2], [1 2])))^2 + ...
             (det(J([2 3], [1 2])))^2 + ...
             (det(J([3 1], [1 2])))^2);
A = int(int(NormJ, xi, -1, 1), eta, -1, 1);
var = vpa(A, 4);

function dn = quad4_deriv (xi, eta)
n = [1.0/4.0 * (1 - xi) * (1 - eta)
     1.0/4.0 * (1 + xi) * (1 - eta)
     1.0/4.0 * (1 + xi) * (1 + eta)
     1.0/4.0 * (1 - xi) * (1 + eta)];
dn = [diff(n, xi), diff(n, eta)];
end

