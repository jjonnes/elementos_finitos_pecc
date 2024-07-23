Exercises 7.2
1. Solve the boundary value problem stated below using the three studied methods of weighted residuals and one basis function
Make a plot with the three approximate functions obtained and a second plot with the three corresponding residual functions.

syms x

%Weighted Residuals
syms c1
N1 = x*(x-1);
y  = @(x) c1*N1;
R  = @(x)diff((y(x)),2) + y(x) - 2*x;
w1 = @(x) 1;
A = int(R(x)*w1(x), 0, 1);
c1 = solve (A, c1);
y  = @(x) c1*N1;
R  = @(x)diff((y(x)),2) + y(x) - 2*x;
y1x = y(x);
R1x = R(x);

%Least squares
syms c1
N1 = x*(x-1);
y  = @(x) c1*N1;
R  = @(x)diff((y(x)),2) + y(x) - 2*x;
w1 = diff(R(x),c1);
A = int(R(x)*w1, 0, 1);
c1 = solve (A, c1);
y  = @(x) c1*N1;
R  = @(x)diff((y(x)),2) + y(x) - 2*x;
y2x = y(x);
R2x = R(x);

%Galerkin
syms c1
N1 = x*(x-1);
y  = @(x) c1*N1;
R  = @(x)diff((y(x)),2) + y(x) - 2*x;
w1 = N1;
A = int(R(x)*w1, 0, 1);
c1 = solve (A, c1);
y  = @(x) c1*N1;
R  = @(x)diff((y(x)),2) + y(x) - 2*x;
y3x = y(x);
R3x = R(x);

X = 0:0.01:1;
y1 = subs (y1x,x, X);
y2 = subs (y2x,x, X);
y3 = subs (y3x,x, X);
R1 = subs (R1x,x, X);
R2 = subs (R2x,x, X);
R3 = subs (R3x,x, X);
plot(X,y1,'-', X,y2,':', X,y3,'-.')
legend ('y1(x)','y2(x)','y3(x)')
plot(X,R1,'-', X,R2,':', X,R3,'-.')
legend ('R1(x)','R2(x)','R3(x)')

2. Solve the problem above using the Galerkin method with two basis functions. You may need symbolic manipulation software inorder to speed up your calculations.

N1 = x*(x-1);
N2 = x*x*(x-1);
y  = @(x) c1*N1+c2*N2;
R  = @(x)diff((y(x)),2) + y(x) - 2*x;
w1 = N1;
w2 = N2;
A1 = int(R(x)*w1, 0, 1);
A2 = int(R(x)*w2, 0, 1);
A = solve ([A1 A2], [c1 c2]);
c1 = A.c1;
c2 = A.c2;
y  = @(x) c1*N1+c2*N2;
R  = @(x)diff((y(x)),2) + y(x) - 2*x;
yx = y(x);
Rx = R(x);

X = 0:0.01:1;
y = subs (yx,x, X);
R = subs (Rx,x, X);
plot(X,y,X,R)
legend ('y(x)','R(x)')

3. Solve the boundary value problem using the Galerkin method and two basis function.
Compare the results with the analytical solution. Explain why the Galerkin method was able to match the analytical result.

syms x
syms c1

N1 = x*(x-1);
N2 = 1-x;
u  = @(x) c1*N1+N2;
R  = @(x) diff((u(x)),2) - diff((u(x))) - 2*x;
w1 = N1;
A = int(R(x)*w1, 0, 1);
c1 = solve (A, c1);
u  = @(x) c1*N1+N2;
R  = @(x) diff((u(x)),2) - diff((u(x))) - 2*x;
ux = u(x);
Rx = R(x);

X = 0:0.01:1;
u = subs (ux,x, X);
R = subs (Rx,x, X);
plot(X,u,X,R)
legend ('y(x)','R(x)')