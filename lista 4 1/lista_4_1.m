%% Exercises 7.1
%1. The nodal displacements at the three nodes of a 1D element are shown in the table. Find the displacement atξ= 0.5.Explain why the result of the interpolation is negative if all nodal values are positive?
n = 101;
xi = linspace (-1, 1, n);      
u = [
0.99,                  % x1 = -1
0.11,                  % x1 = 1
0.05];                 % x1 = 0
for num = 1:n
    N = shape_form_1d(xi(num));
    U(num,1) = N * u;
end
figure
plot(xi,U,[-1,1,0],u,'*')
grid;
xlabel ("ξ");
ylabel ("u");

function N = shape_form_1d(xi)
N = [
1/2*xi^2 - 1/2*xi;     % N1
1/2*xi^2 + 1/2*xi;     % N2
1 - xi^2               % N3
]';                   
end
%% Exercises 7.1
%2.The values of a nodal variableuat the five nodes of a 1D element are shown in the table below. Plot theufield along theelement by using shape functions.
n = 101;
xi = linspace (-1, 1, n);      
u = [0.00; 0.22; -0.15; 0.12; 0.15];               
for num = 1:n
    N = shape_form_1d(xi(num));
    U(num,1) = N * u;
end
figure
plot(xi,U,[-1,1,-1/2,1/2,0],u,'*')
grid;
xlabel ("ξ");
ylabel ("u");

function N = shape_form_1d(xi)
N = [
1/6*(    xi -  xi^2 - 4*xi^3 + 4*xi^4);   % N1    
1/6*(-   xi -  xi^2 + 4*xi^3 + 4*xi^4);   % N2
1/3*(-4*xi + 8*xi^2 + 4*xi^3 - 8*xi^4);   % N3
1/3*( 4*xi + 8*xi^2 - 4*xi^3 - 8*xi^4);   % N4
(1         - 5*xi^2          + 4*xi^4)    % N5
]';           
end
%% Exercises 7.1
%3. Determine the shape functions for a 1D isoparametric element with four and five nodes by calculating the polynomialcoefficients using a system of equations.  Later, compare the results with the shape functions obtained using Lagrangepolynomials. Plot the obtained shape functions.

5 Nós
n = 101;
xi = linspace (-1, 1, n);      
u1 = [1.00 0.00 0.00 0.00 0.00]';
u2 = [0.00 1.00 0.00 0.00 0.00]';               
u3 = [0.00 0.00 1.00 0.00 0.00]';               
u4 = [0.00 0.00 0.00 1.00 0.00]';               
u5 = [0.00 0.00 0.00 0.00 1.00]';               

for i = 1:n
    N5 = shape_form_1d_5n(xi(i));
    U1(i, 1) =  N5 * u1;

    N5 = shape_form_1d_5n(xi(i));
    U2(i, 1) =  N5 * u2;

    N5 = shape_form_1d_5n(xi(i));
    U3(i, 1) =  N5 * u3;

    N5 = shape_form_1d_5n(xi(i));
    U4(i, 1) =  N5 * u4;

    N5 = shape_form_1d_5n(xi(i));
    U5(i, 1) =  N5 * u5;
end

figure
plot(xi,U1,[-1,1,-1/2,1/2,0],u1,'*')
grid;
xlabel ("ξ");
hold on
plot(xi,U2,[-1,1,-1/2,1/2,0],u2,'*')
hold on
plot(xi,U3,[-1,1,-1/2,1/2,0],u3,'*')
hold on
plot(xi,U4,[-1,1,-1/2,1/2,0],u4,'*')
hold on
plot(xi,U5,[-1,1,-1/2,1/2,0],u5,'*')
hold on
clear all

4 Nós
n = 101;
xi = linspace (-1, 1, n);      
u1 = [1.00 0.00 0.00 0.00]';
u2 = [0.00 1.00 0.00 0.00]';               
u3 = [0.00 0.00 1.00 0.00]';               
u4 = [0.00 0.00 0.00 1.00]';               

for i = 1:n
    N4 = shape_form_1d_4n(xi(i));
    U1(i, 1) =  N4 * u1;

    N4 = shape_form_1d_4n(xi(i));
    U2(i, 1) =  N4 * u2;

    N4 = shape_form_1d_4n(xi(i));
    U3(i, 1) =  N4 * u3;

    N4 = shape_form_1d_4n(xi(i));
    U4(i, 1) =  N4 * u4;
end

figure
plot(xi,U1,[-1,1,-1/3,1/3],u1,'*')
grid;
xlabel ("ξ");
hold on
plot(xi,U2,[-1,1,-1/3,1/3],u2,'*')
hold on
plot(xi,U3,[-1,1,-1/3,1/3],u3,'*')
hold on
plot(xi,U4,[-1,1,-1/3,1/3],u4,'*')
hold on


function N5 = shape_form_1d_5n(xi)
N5 = [
1/6*(    xi -  xi^2 - 4*xi^3 + 4*xi^4);   % N1    
1/6*(-   xi -  xi^2 + 4*xi^3 + 4*xi^4);   % N2
1/3*(-4*xi + 8*xi^2 + 4*xi^3 - 8*xi^4);   % N3
1/3*( 4*xi + 8*xi^2 - 4*xi^3 - 8*xi^4);   % N4
(1         - 5*xi^2          + 4*xi^4)    % N5
]';           
end

function N4 = shape_form_1d_4n(xi)
N4 = [
1/16*(-1 +    xi + 9*xi^2 -  9*xi^3);    % N1    
1/16*(-1 -    xi + 9*xi^2 +  9*xi^3);    % N2
1/16*( 9 - 27*xi - 9*xi^2 + 27*xi^3);    % N3
1/16*( 9 + 27*xi - 9*xi^2 - 27*xi^3)     % N4
]';
         
end
%% Exercises 7.1
%4. Using computer software elaborate a routine that computes the shape function coefficients of an-node 1D element, with n= 2,3,4,5.
clc
clear all

syms xi
n = input('Número de nós: ');
x = linspace(-1,1,n);

if rem(n,2) ~= 0
impar = linspace (1,n,(ceil(n/2)));
par = linspace (n-1,2,(floor(n/2)));
else
impar = linspace (1,n-1,((n/2)));
par = linspace (n,2,((n/2)));
end
k = [impar par];

for num = 1:n
    N(k(num),1) = shapefun (x, num, xi, n);
end

simplifyFraction(simplify(N))

function ni = shapefun (x, i, xi, n)    
for j = 1:n
    if i~=j
       a(j,1) = (xi-x(j)) / (x(i)-x(j));
    for j = 1:n
       if i==j
          a(j,1) = 1;
       end
    end
    end
end
ni = prod(a);
end
