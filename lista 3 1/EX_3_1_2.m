%%Exercises 3.1 - 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Determine the nodal displacments for the second truss shown
%%at rigth. Use P = 10 kN and EA = 2000 kN. The lengths for
%%elements 1 and 5 are 1 m and 0.8 m, respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
format short
format compact
close all


%% Dados da estrutura
%MN/m�
E = 2000;                       
%m�
A = 1;                        


%% Geometria
% N�s [NumNo, Cx, Cy] [m]
Nos = [1 0 0; 2 1 0; 3 1 0.8; 4 0 0.8];


%% Elementos
% Barra [NumElem, No'i', No'j']
Barra = [1 1 2; 2 4 3; 3 1 3; 4 4 2; 5 2 3];


%% Restri�oes de apoio
% Restri��es de apoio [NumNo, Rx, Ry]
Restr = [1 1 1; 4 1 1];


%% Carregamento
% Carregametno [No, Px, Py]
P = [3 0 -10];


%% Montagem inicial de matrizes e vetores
% Grau de liberdade
GL = size(Nos, 1) * 2;
% N�mero de elementos 
NumElem = size(Barra, 1);
% Matriz K (elemento)
Ke = zeros(2, 2, NumElem);
% Matriz K (Global)
K = zeros(GL, GL);
% Vetor U
U = zeros(GL, 1);
% Vetor F
F = zeros(GL, 1);
% For�a de rea��o
Fr = zeros(GL, 1);
% Comprimento do elemnto
L = zeros(NumElem, 1);
% Coordenadas por N�
Coor = zeros(NumElem, 4);
% Seno e coseno
CS = zeros(NumElem,2);


%% Matriz de rigidez K (Global)
for elem = 1:NumElem
    % N�s 'i' 'j'
    NoI = Barra(elem, 2);
    NoJ = Barra(elem, 3);
    % Determina��o L
    L(elem) = norm(Nos(NoJ, [2, 3]) - Nos(NoI, [2, 3]));
    % Determina��o das coordenadas nodais 
    Coor(elem, :) = [NoI * 2 - [1, 0], NoJ * 2 - [1, 0]];
    % Matriz de transforma��o T
    CS(elem, 1) = (Nos(NoJ, 2) - Nos(NoI, 2)) /L(elem); 
    CS(elem, 2) = (Nos(NoJ, 3) - Nos(NoI, 3)) /L(elem);
    C = CS(elem, 1); S = CS(elem, 2);
    T = [C S 0 0; 0 0 C S];
    % Determina��o Ke
    Ke(:, :, elem) = (E * A / L(elem)) * [1  -1; -1 1];
    % Determina��o K
    K(Coor(elem, :) , Coor(elem, :)) = K(Coor(elem, :) , Coor(elem, :)) + T' *  Ke(:, :, elem) * T;
    KT(:, :, elem) = T' * [1 -1;-1 1] * T;
end


%% Vetor de for�as
for load = 1: size(P, 1)
    % N�s
    NoP = P(load, 1);
    % Determina��o das cordenadas 
    CoorP = (NoP * 2 - [1, 0])';
    % For�a aplicada
    Fa = P(load, [2, 3])';
    % Vetor de for�as
    F(CoorP) = Fa;   
    
end


%% Restri��es nos apoios
% Determina��o do n�mero de restri��es
NumGLR = sum(sum(Restr(:, [2, 3])));
% Vetor deslocamentos restringidos 
GLR = zeros(NumGLR, 1);


%% Cordenadas com resti�oes
i = 1;
for apoio = 1:(size(Restr, 1))
    % N� de apoio
    NoApoio = Restr(apoio, 1);
    % Restri��o em x
    if Restr(apoio, 2) == 1
        GLR(i, 1) = NoApoio * 2 - 1;
        i = i + 1;
    end
    % Restri��o em y
    if Restr(apoio, 3) == 1
        GLR(i, 1) = NoApoio * 2;
        i = i + 1;
    end
end


%% Cordenadas sem restri��o
% Grau de liberdade sem restri��o
GLSR = setxor((1:GL)', GLR);


%% Deslocamentos sem resti��o
% Deslocamentos
U(GLSR) = K(GLSR, GLSR) \ F(GLSR);


%% For�as de rea��o
% For�a de rea��o
Fr(GLR) = K(GLR, :) * U;


%% Tensoes  e For�as Normais
FNormal = zeros(NumElem ,1);
SigmaX = zeros(NumElem, 1);
for elem = 1:NumElem
    C = CS(elem, 1); S = CS(elem, 2);
    T = [C S 0 0; 0 0 C S];
    UAxial(:, :, elem) = T * U(Coor(elem, :));
    DeltaL(elem) = UAxial(2, 1, elem) - UAxial(1, 1, elem);
    SigmaX(elem, 1) = E / L(elem) * DeltaL(elem); 
    FNormal(elem) = (DeltaL(elem) * A * E) / L(elem);
end
SigmaX;
FNormal;