%% Exercises 4.1
%1. Find the Jacobian determinants for the transformation from Cartesian coordinates to polar coordinates
clear all
clc

syms x y r theta phi rho

x = r * cos(theta);
y = r * sin(theta);

jacobian1 = [diff(x, r) diff(x, theta)
            diff(y, r) diff(y, theta)];

J1 = simplify(det(jacobian1));

clear x y 

x = rho * sin(phi) * cos(theta);
y = rho * sin(phi) * sin(theta);
z = rho * cos(phi);

jacobian2 = [diff(x, rho) diff(x, theta) diff(x, phi)
             diff(y, rho) diff(y, theta) diff(y, phi)
             diff(z, rho) diff(z, theta) diff(z, phi)];

J2 = simplify(det(jacobian2));
%% Exercises 4.1
%2. Calculate the integral below and show the Jacobian matrix used.I=∫+∞−∞∫+∞−∞e−x2−y2dxdy
clear all 
clc

syms x y r theta

x = r * cos(theta);
y = r * sin(theta);

jacobian = [diff(x, r) diff(x, theta)
            diff(y, r) diff(y, theta)];

J = simplify(det(jacobian));

eqxy = exp(-x^2-y^2);
eqrtheta = eqxy * J;

A = int(int(eqrtheta, r, 0, inf), theta, 0, 2*pi);
%% Exercises 4.1
%3. Nodal coordinates of a 1D element are:x1= 2.0,x2= 10.0andx3= 9.0. Plot the JacobianJon the natural coordinatesspace and explain the behavior between nodes 2 and 3.
clear all
clc

n = 101;
xi = linspace (-1, 1, n);

C = [2; 10; 9];
for i = 1:n
    J(i, 1) = C'* lin_deriv (xi(i));
end
plot(xi, J)
ylabel ("J")
xlabel ("ξ")
grid on

function dn = lin_deriv (xi)
dn = [xi - 1/2
      xi + 1/2
      -2*xi];
end
%% Exercises 4.1
%5. Write a routine that given the coordinates of a quadrilateral element with four nodes calculates the Jacobian at givencoordinates(ξ,η).

clear all
clc

C = input('Coordenadas (x y):');
xi = input ('xi:');
eta = input ('eta:');
jacobian = C' * quad4_derivs (xi, eta);
J = det(jacobian);

function dn = quad4_derivs(xi, eta)
dn = [0.25*(-1.0+eta)   0.25*(-1.0+xi)
      0.25*(+1.0-eta)   0.25*(-1.0-xi)
      0.25*(+1.0+eta)   0.25*(+1.0+xi)
      0.25*(-1.0-eta)   0.25*(+1.0-xi)];
end

