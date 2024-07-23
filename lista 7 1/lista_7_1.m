1. Find an approximate solution for the BVP using weighted integrals with an approximate solution of the form u(x) =c1x2+c2x+c3 and a weight function w(x) = 1.

syms c1 c2 c3
syms x

u  = @(x) c1*x*x+c2*x+c3;
w1 = @(x) 1;
w2 = @(x) x;
R  = @(x) diff(u(x),2) + u(x);
u0 = u(0) == 1;
u1 = u(1) == 0;
A1  = int(R(x)*w1(x), 0, 1) == 1;
A2  = int(R(x)*w2(x), 0, 1) == 1;

eqs = [subs(A1, c3, 1), ...
       subs(A2, c3, 1)];
A = solve (eqs, [c1 c2]);
C1 = A.c1;
C2 = A.c2;
C3 = 1;
u  = @(x) C1*x*x+C2*x+C3;
R  = @(x) diff(u(x),2) + u(x);
ux = u(x);
Rx = R(x);


2. Find an approximate solution for the BVP using weighted integrals with an approximate solution of the form u(x) =c1x3+c2x2+c3x+c4 and weight function sw1(x) = 1 and w2(x) =x.

syms c1 c2 c3 c4
syms x
u  = @(x) c1*x^3+c2*x^2+c3*x+c4;
w1 = @(x) 1;
w2 = @(x) x;
R  = @(x) diff(diff(u(x))*x) + u(x);
u0 = u(0) == 1;
ux = diff(u(x))*x == 0;
ux1 = subs (ux, x, 1);
A1  = int(R(x)*w1(x), 0, 1) == 0;
A2  = int(R(x)*w2(x), 0, 1) == 0;

eqs = [subs(A1, c4, 1), ...
       subs(A2, c4, 1)];
A = solve ([eqs ux1], [c1 c2 c3]);
C1 = A.c1;
C2 = A.c2;
C3 = A.c3;
C4 = 1;
u  = @(x) C1*x^3+C2*x^2+C3*x+C4;
R  = @(x) diff(diff(u(x))*x) + u(x);
ux = u(x);
Rx = R(x);




