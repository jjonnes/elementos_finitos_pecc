%% Exercises 7.2
%2. List all terms needed to generate the shape functions for a serendipity quadrilateral element with 12 nodes. Use computersoftware to compute all polynomials coefficients.
clear all 
clc

syms xi eta 
syms a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 
syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 

Cxy = a1 + a2*xi + a3*eta + a4*xi*eta + a5*xi^2 + ...
    a6*eta^2 + a7*xi^2*eta + a8*xi*eta^2 + a9*xi^3 + ... 
    a10*eta^3 + a11*xi^3*eta + a12*xi*eta^3;

x = linspace (-1, 1, 4);

eqs =[subs(Cxy,[xi eta],[x(1) x(1)]) == c1, ...
    subs(Cxy,[xi eta],[x(4) x(1)]) == c2,...
    subs(Cxy,[xi eta],[x(4) x(4)]) == c3,...
    subs(Cxy,[xi eta],[x(1) x(4)]) == c4,...
    subs(Cxy,[xi eta],[x(2) x(1)]) == c5,...
    subs(Cxy,[xi eta],[x(3) x(1)]) == c6,...
    subs(Cxy,[xi eta],[x(4) x(2)]) == c7,...
    subs(Cxy,[xi eta],[x(4) x(3)]) == c8,...
    subs(Cxy,[xi eta],[x(3) x(4)]) == c9,...
    subs(Cxy,[xi eta],[x(2) x(4)]) == c10,...
    subs(Cxy,[xi eta],[x(1) x(3)]) == c11,...
    subs(Cxy,[xi eta],[x(1) x(2)]) == c12];

var = [a1, a2,a3, a4, a5, a6, a7, a8, a9, a10, a11, a12];
Cvar = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12];

A = solve(eqs,var);
a1 = A.a1;
a2 = A.a2;
a3 = A.a3;
a4 = A.a4;
a5 = A.a5;
a6 = A.a6;
a7 = A.a7;
a8 = A.a8;
a9 = A.a9;
a10 = A.a10;
a11 = A.a11;
a12 = A.a12;

Cxy = a1 + a2*xi + a3*eta + a4*xi*eta + a5*xi^2 + ...
    a6*eta^2 + a7*xi^2*eta + a8*xi*eta^2 + a9*xi^3 + ... 
    a10*eta^3 + a11*xi^3*eta + a12*xi*eta^3;

[N,Ci] = coeffs(Cxy,Cvar);
N1  = simplify(N( 1));
N2  = simplify(N( 2));
N3  = simplify(N( 3));
N4  = simplify(N( 4));
N5  = simplify(N( 5));
N6  = simplify(N( 6));
N7  = simplify(N( 7));
N8  = simplify(N( 8));
N9  = simplify(N( 9));
N10 = simplify(N(10));
N11 = simplify(N(11));
N12 = simplify(N(12));

sum (N1+N2+N3+N4+N5+N6+N7+N8+N9+N10+N11+N12);

fsurf(N1,[-1 1 -1 1])
title 'N1'
xlabel ("ξ");
ylabel ("η");
fsurf(N2,[-1 1 -1 1])
title 'N2'
xlabel ("ξ");
ylabel ("η");
fsurf(N3,[-1 1 -1 1])
title 'N3'
xlabel ("ξ");
ylabel ("η");
fsurf(N4,[-1 1 -1 1])
title 'N4'
xlabel ("ξ");
ylabel ("η");
fsurf(N5,[-1 1 -1 1])
title 'N5'
xlabel ("ξ");
ylabel ("η");
fsurf(N6,[-1 1 -1 1])
title 'N6'
xlabel ("ξ");
ylabel ("η");
fsurf(N7,[-1 1 -1 1])
title 'N7'
xlabel ("ξ");
ylabel ("η");
fsurf(N8,[-1 1 -1 1])
title 'N8'
xlabel ("ξ");
ylabel ("η");
fsurf(N9,[-1 1 -1 1])
title 'N9'
xlabel ("ξ");
ylabel ("η");
fsurf(N10,[-1 1 -1 1])
title 'N10'
xlabel ("ξ");
ylabel ("η");
fsurf(N11,[-1 1 -1 1])
title 'N11'
xlabel ("ξ");
ylabel ("η");
fsurf(N12,[-1 1 -1 1])
title 'N12'
xlabel ("ξ");
ylabel ("η");
%% Exercises 7.2
%3. Find the shape functions of nodes 1, 4 and 10 of a 10-node 2D Lagrangian triangular element.
clear all 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       1  3                  %
%       |  |\                 %
%       |  | \                %
%       |  8  7               %
%       |  |   \              %
%      eta |    \             %
%       |  9 10  6            %
%       |  |      \           %
%       |  |       \          %
%       0  1--4--5--2         %
%       |                     %
%       +--0----xi---1-->     % 
%                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


syms xi eta 
syms a1 a2 a3 a4 a5 a6 a7 a8 a9 a10
syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10

Cxy = a1 + a2*xi + a3*eta + ...
    a4*xi^2 + a5*eta^2 + a6*xi*eta + ...
    a7*xi^3 + a8*eta^3 + a9*xi^2*eta + a10*xi*eta^2;

x = linspace (0, 1, 4);

eqs =[subs(Cxy,[xi eta],[x(1) x(1)]) == c1, ...
    subs(Cxy,[xi eta],[x(4) x(1)]) == c2,...
    subs(Cxy,[xi eta],[x(1) x(4)]) == c3,...
    subs(Cxy,[xi eta],[x(2) x(1)]) == c4,...
    subs(Cxy,[xi eta],[x(3) x(1)]) == c5,...
    subs(Cxy,[xi eta],[x(3) x(2)]) == c6,...
    subs(Cxy,[xi eta],[x(2) x(3)]) == c7,...
    subs(Cxy,[xi eta],[x(1) x(3)]) == c8,...
    subs(Cxy,[xi eta],[x(1) x(2)]) == c9,...
    subs(Cxy,[xi eta],[x(2) x(2)]) == c10];

var = [a1, a2,a3, a4, a5, a6, a7, a8, a9, a10];
Cvar = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10];

A = solve(eqs,var);
a1 = A.a1;
a2 = A.a2;
a3 = A.a3;
a4 = A.a4;
a5 = A.a5;
a6 = A.a6;
a7 = A.a7;
a8 = A.a8;
a9 = A.a9;
a10 = A.a10;

Cxy = a1 + a2*xi + a3*eta + ...
    a4*xi^2 + a5*eta^2 + a6*xi*eta + ...
    a7*xi^3 + a8*eta^3 + a9*xi^2*eta + a10*xi*eta^2;

[N,Ci] = coeffs(Cxy,Cvar);
N1  = simplify(N( 1));
N2  = simplify(N( 2));
N3  = simplify(N( 3));
N4  = simplify(N( 4));
N5  = simplify(N( 5));
N6  = simplify(N( 6));
N7  = simplify(N( 7));
N8  = simplify(N( 8));
N9  = simplify(N( 9));
N10 = simplify(N(10));

sum (N(1)+N(2)+N(3)+N(4)+N(5)+N(6)+N(7)+N(8)+N(9)+N(10));

fsurf(N1,[0 1 0 1])
title 'N1'
xlabel ("ξ");
ylabel ("η");
fsurf(N2,[0 1 0 1])
title 'N2'
xlabel ("ξ");
ylabel ("η");
fsurf(N3,[0 1 0 1])
title 'N3'
xlabel ("ξ");
ylabel ("η");
fsurf(N4,[0 1 0 1])
title 'N4'
xlabel ("ξ");
ylabel ("η");
fsurf(N5,[0 1 0 1])
title 'N5'
xlabel ("ξ");
ylabel ("η");
fsurf(N6,[0 1 0 1])
title 'N6'
xlabel ("ξ");
ylabel ("η");
fsurf(N7,[0 1 0 1])
title 'N7'
xlabel ("ξ");
ylabel ("η");
fsurf(N8,[0 1 0 1])
title 'N8'
xlabel ("ξ");
ylabel ("η");
fsurf(N9,[0 1 0 1])
title 'N9'
xlabel ("ξ");
ylabel ("η");
fsurf(N10,[0 1 0 1])
title 'N10'
xlabel ("ξ");
ylabel ("η");
%% Exercises 7.2
%4. Compute all polynomial coefficients of the shape functions of a 10-node tetrahedron element.
clear all 
clc

syms xi eta zeta
syms a1 a2 a3 a4 a5 a6 a7 a8 a9 a10
syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10

Cxyz = a1 + a2*xi + a3*eta + a4*zeta + ...
    a5*xi^2 + a6*eta^2 + a7*zeta^2 + ... 
    a8*xi*zeta + a9*eta*zeta + a10*xi*eta;

x = linspace (0, 1, 3);

eqs =[subs(Cxyz,[xi eta zeta],[x(1) x(1) x(1)]) == c1, ...
    subs(Cxyz,[xi eta zeta],[x(3) x(1) x(1)]) == c2,...
    subs(Cxyz,[xi eta zeta],[x(1) x(3) x(1)]) == c3,...
    subs(Cxyz,[xi eta zeta],[x(1) x(1) x(3)]) == c4,...
    subs(Cxyz,[xi eta zeta],[x(2) x(1) x(1)]) == c5,...
    subs(Cxyz,[xi eta zeta],[x(2) x(2) x(1)]) == c6,...
    subs(Cxyz,[xi eta zeta],[x(1) x(2) x(1)]) == c7,...
    subs(Cxyz,[xi eta zeta],[x(1) x(1) x(2)]) == c8,...
    subs(Cxyz,[xi eta zeta],[x(2) x(1) x(2)]) == c9,...
    subs(Cxyz,[xi eta zeta],[x(1) x(2) x(2)]) == c10];

var = [a1, a2,a3, a4, a5, a6, a7, a8, a9, a10];
Cvar = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10];

A = solve(eqs,var);
a1 = A.a1;
a2 = A.a2;
a3 = A.a3;
a4 = A.a4;
a5 = A.a5;
a6 = A.a6;
a7 = A.a7;
a8 = A.a8;
a9 = A.a9;
a10 = A.a10;

Cxyz = a1 + a2*xi + a3*eta + a4*zeta + ...
    a5*xi^2 + a6*eta^2 + a7*zeta^2 + ... 
    a8*xi*eta + a9*eta*zeta + a10*xi*zeta;

[N,Ci] = coeffs(Cxyz,Cvar);
N1  = simplify(N( 1));
N2  = simplify(N( 2));
N3  = simplify(N( 3));
N4  = simplify(N( 4));
N5  = simplify(N( 5));
N6  = simplify(N( 6));
N7  = simplify(N( 7));
N8  = simplify(N( 8));
N9  = simplify(N( 9));
N10 = simplify(N(10));

sum (N(1)+N(2)+N(3)+N(4)+N(5)+N(6)+N(7)+N(8)+N(9)+N(10))
%% Exercises 7.2
%5. (Opt.) Find all shape functions for the transition element shown below.
%6. (Opt.) Using computer software plot the shape functions of nodes 1 and 5 for the element in the last exercise.
clear all 
clc

syms xi eta a b x

% ELEMENTO DE TRANSIÇÃO QUADRILATERAL COM 5 NÓS
NC =  1/4 * (1 + xi*a) * (1 + eta*b);    % i = 1...4 (linear)
NMP = 1/2 * (1 + xi*a) * (1 - eta^2);    % i = 6, 8  (quadratico)
NMI = 1/2 * (1 + xi^2) * (1 - eta*b);    % i = 5, 7  (quadratico)

N1 = subs (NC,  [a b], [-1 -1]);
N2 = subs (NC,  [a b], [ 1 -1]);
N3 = subs (NC,  [a b], [ 1  1]);
N4 = subs (NC,  [a b], [-1  1]);
N5 = subs (NMP, [a b], [ 1  0]);
k = sum(N1+N2+N3+N4+N5);


fsurf(N1,[-1 1 -1 1])
title 'N1'
xlabel ("ξ");
ylabel ("η");
fsurf(N2,[-1 1 -1 1])
title 'N2'
xlabel ("ξ");
ylabel ("η");
fsurf(N3,[-1 1 -1 1])
title 'N3'
xlabel ("ξ");
ylabel ("η");
fsurf(N4,[-1 1 -1 1])
title 'N4'
xlabel ("ξ");
ylabel ("η");
fsurf(N5,[-1 1 -1 1])
title 'N5'
xlabel ("ξ");
ylabel ("η");
