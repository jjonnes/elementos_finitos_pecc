%% Exercises 3.2 - 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Compute the nodal displacements and reaction forces in the cantilever
%%beam. Use two finite elements and assume EI = 100 kNm2,
%%L = 1 m, M = 20 KNm and q = 120 kN/m. Draw the diagrams for
%%shear forces, bending moment and deflection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
format short
format compact

%% Dados da estrutura
%MN/m²
E = 100;                       
%m²
I = 1;                        

%% Geometria
% Nós [NumNo, Cx, Cy] [m]
Nos = [1 0 0; 2 1 0; 3 2 0];

%% Elementos
% Barra [NumElem, No'i', No'j']
Barra = [1 1 2; 2 2 3];

%% Restriçoes de apoio
% Restrições de apoio [NumNo, R, M]
Restr = [1 1 1];

%% Carregamento
% Carregametno [No, P, M]
P = [2 0 20];

%% Montagem inicial de matrizes e vetores
% Grau de liberdade
GL = size(Nos, 1) * 2;
% Número de elementos 
NumElem = size(Barra, 1);
% Matriz K (elemento)
Ke = zeros(4, 4, NumElem);
% Matriz K (Global)
K = zeros(GL, GL);
% Vetor U
U = zeros(GL, 1);
% Vetor F
F = zeros(GL, 1);
% Força de reação
Fr = zeros(GL, 1);
% Comprimento do elemnto
L = zeros(NumElem, 1);
% Coordenadas por Nó
Coor = zeros(NumElem, 4);
% Seno e coseno
CS = zeros(NumElem,2);

%% Matriz de rigidez K (Global)
for elem = 1:NumElem
    % Nós 'i' 'j'
    NoI = Barra(elem, 2);
    NoJ = Barra(elem, 3);
    % Determinação L
    L(elem) = norm(Nos(NoJ, [2, 3]) - Nos(NoI, [2, 3]));
    % Determinação das coordenadas nodais 
    Coor(elem, :) = [NoI * 2 - [1, 0], NoJ * 2 - [1, 0]];
    % Determinação Ke
    Ke(:, :, elem) = (E * I / (L(elem))^3) * [12  6*L(elem) -12 6*L(elem);
        6*L(elem) 4*(L(elem))^2 -6*L(elem) 2*(L(elem))^2; -12  -6*L(elem) 12 -6*L(elem); 
        6*L(elem) 2*(L(elem))^2 -6*L(elem) 4*(L(elem))^2];
    % Determinação K
    K(Coor(elem, :) , Coor(elem, :)) = K(Coor(elem, :) , Coor(elem, :)) +  Ke(:, :, elem);
end

%% Vetor de forças
for load = 1: size(P, 1)
    % Nós
    NoP = P(load, 1);
    % Determinação das cordenadas 
    CoorP = (NoP * 2 - [1, 0])';
    % Força aplicada
    Fa = P(load, [2, 3])';
    % Vetor de forças
    F(CoorP) = Fa;   
end

%% Restrições nos apoios
% Determinação do número de restrições
NumGLR = sum(sum(Restr(:, [2, 3])));
% Vetor deslocamentos restringidos 
GLR = zeros(NumGLR, 1);

%% Cordenadas com restiçoes
i = 1;
for apoio = 1:(size(Restr, 1))
    % Nó de apoio
    NoApoio = Restr(apoio, 1);
    % Restrição em x
    if Restr(apoio, 2) == 1
        GLR(i, 1) = NoApoio * 2 - 1;
        i = i + 1;
    end
    % Restrição em y
    if Restr(apoio, 3) == 1
        GLR(i, 1) = NoApoio * 2;
        i = i + 1;
    end
end

%% Cordenadas sem restrição
% Grau de liberdade sem restrição
GLSR = setxor((1:GL)', GLR);

%% Deslocamentos sem restição
% Deslocamentos
U(GLSR) = K(GLSR, GLSR) \ (F(GLSR) + [-42 -6 -18 4]');

%% Forças de reação
% Força de reação
Fr(GLR) = K(GLR, :) * U;

%% Diagramas
x = 1;

l = 1;
n1 = (1 - 3*x^2/l^2 + 2*x^3/l^3);
n2 = (x - 2*x^2/l   +   x^3/l^2);
n3 = (  3*x^2/l^2 - 2*x^3/l^3);
n4 = (   x^3/l^2 -   x^2/l );

n11  = (-6*x + 6*x^2);
n12  = (1 - 4*x + 3*x^2);
n13  = (6*x - 6*x^2);
n14  = (-2*x + 3*x^2);

n21  = (-6 + 12*x);
n22  = (-4 + 6*x);
n23  = (6 - 12*x);
n24  = (-2 + 6*x);

n31  = 12;
n32  = 6;
n33  = -12;
n34  = 6;

a1 = [n1 n2 n3 n4] * U([Coor(2,:)]);
a2 = [n11 n12 n13 n14] * U([Coor(2,:)]);
a3 = [n21 n22 n23 n24] * U([Coor(2,:)]);
a4 = [n31 n32 n33 n34] * U([Coor(2,:)]);


%% Exercises 3.2 - 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Compute the displacements and reaction forces in the one story
%%portal frame. Consider P = 10 kN, EI = 10 kNm2, L = 1 m and
%%EA = 100 kN. Draw the diagrams for bending moment and shear
%%forces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
format short
format compact

%% Dados da estrutura
%MN/m²
E = 1;                       
%m²
I = 10;   

A = 100;

%% Geometria
% Nós [NumNo, Cx, Cy] [m]
Nos = [1 0 1; 2 1 1; 3 0 0; 4 1 0];


%% Elementos
% Barra [NumElem, No'i', No'j']
Barra = [1 1 3; 2 1 2; 3 2 4];


%% Restriçoes de apoio
% Restrições de apoio [NumNo, Rx, Ry, M]
Restr = [3 1 1 1; 4 1 1 1];

%% Carregamento
% Carregametno [No, Px, Py, M]
P = [1 10 0 0];

%% Montagem inicial de matrizes e vetores
% Grau de liberdade
GL = size(Nos, 1) * 3;
% Número de elementos 
NumElem = size(Barra, 1);
% Matriz K (elemento)
Ke = zeros(6, 6, NumElem);
KeBeam = zeros(6, 6, NumElem);
KeBarra = zeros(6, 6, NumElem);
% Matriz K (Global)
K = zeros(GL, GL);
% Vetor U
U = zeros(GL, 1);
% Vetor F
F = zeros(GL, 1);
% Força de reação
Fr = zeros(GL, 1);
% Comprimento do elemnto
L = zeros(NumElem, 1);
% Coordenadas por Nó
Coor = zeros(NumElem, 6);
% Seno e coseno
CS = zeros(NumElem,2);

%% Matriz de rigidez K (Global)
for elem = 1:NumElem
    % Nós 'i' 'j'
    NoI = Barra(elem, 2);
    NoJ = Barra(elem, 3);
    % Determinação L
    L(elem) = norm(Nos(NoJ, [2, 3]) - Nos(NoI, [2, 3]));
    % Determinação das coordenadas nodais 
    Coor(elem, :) = [NoI * 3 - [2, 1, 0], NoJ * 3 - [2, 1, 0]];
    % Matriz de transformação T
    CS(elem, 1) = (Nos(NoJ, 2) - Nos(NoI, 2)) /L(elem); 
    CS(elem, 2) = (Nos(NoJ, 3) - Nos(NoI, 3)) /L(elem);
    C = CS(elem, 1); S = CS(elem, 2);
    T = [C S 0 0 0 0; -S C 0 0 0 0; 0 0 1 0 0 0; 0 0 0 C S 0; 0 0 0 -S C 0; 0 0 0 0 0 1];
    % Determinação Ke
    KeBeam([2 3 5 6], [2 3 5 6], elem) = (E * I / (L(elem))^3) * [12  6*L(elem) -12 6*L(elem); 6*L(elem) 4*(L(elem))^2 -6*L(elem) 2*(L(elem))^2;
        -12  -6*L(elem) 12 -6*L(elem); 6*L(elem) 2*(L(elem))^2 -6*L(elem) 4*(L(elem))^2];
    KeBarra([1 4], [1 4], elem) = (E * A / L(elem)) * [1  -1; -1 1];
    Ke = KeBeam + KeBarra;
    % Determinação K
    K(Coor(elem, :) , Coor(elem, :)) = K(Coor(elem, :) , Coor(elem, :)) + T' * Ke(:, :, elem) * T;
    KT(:, :, elem) = T' * (Ke(:, :, elem)) * T;
end

%% Vetor de forças
for load = 1: size(P, 1)
    % Nós
    NoP = P(load, 1);
    % Determinação das cordenadas 
    CoorP = NoP * 3 - [2, 1, 0]';
    % Força aplicada
    Fa = P(load, [2, 3, 4])';
    % Vetor de forças
    F(CoorP) = Fa;   
end

%% Restrições nos apoios
% Determinação do número de restrições
NumGLR = sum(sum(Restr(:, [2, 3, 4])));
% Vetor deslocamentos restringidos 
GLR = zeros(NumGLR, 1);

%% Coordenadas com restiçoes
i = 1;
for apoio = 1:(size(Restr, 1))
    % Nó de apoio
    NoApoio = Restr(apoio, 1);
    % Restrição em x
    if Restr(apoio, 2) == 1
        GLR(i, 1) = NoApoio * 3 - 2;
        i = i + 1;
    end
    % Restrição em y
    if Restr(apoio, 3) == 1
        GLR(i, 1) = NoApoio * 3 - 1;
        i = i + 1;
    end
     % Restrição em M
    if Restr(apoio, 4) == 1
        GLR(i, 1) = NoApoio * 3;
        i = i + 1;
    end
end

%% Cordenadas sem restrição
% Grau de liberdade sem restrição
GLSR = setxor((1:GL)', GLR);

%% Deslocamentos sem restição
% Deslocamentos
U(GLSR) = K(GLSR, GLSR) \ F(GLSR);

%% Forças de reação
% Força de reação
Fr(GLR) = K(GLR, :) * U;


%% Diagramas
x = 0;

l = 1;
n1 = (1 - 3*x^2/l^2 + 2*x^3/l^3);
n2 = (x - 2*x^2/l   +   x^3/l^2);
n3 = (  3*x^2/l^2 - 2*x^3/l^3);
n4 = (   x^3/l^2 -   x^2/l );

n11  = (-6*x + 6*x^2);
n12  = (1 - 4*x + 3*x^2);
n13  = (6*x - 6*x^2);
n14  = (-2*x + 3*x^2);

n21  = (-6 + 12*x);
n22  = (-4 + 6*x);
n23  = (6 - 12*x);
n24  = (-2 + 6*x);

n31  = 12;
n32  = 6;
n33  = -12;
n34  = 6;

a1 = [0 n1 n2 0 n3 n4] * U([Coor(3,:)]); 
a2 = [0 n11 n12 0 n13 n14] * U([Coor(3,:)]); 
a3 = [0 n21 n22 0 n23 n24] * U([Coor(3,:)]);
a4 = [0 n31 n32 0 n33 n34] * U([Coor(3,:)]);
%%

