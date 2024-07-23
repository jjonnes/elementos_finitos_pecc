%%Exercises 3.1 - 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
format short
format compact



%% Dados da estrutura
%kN/m
Rigidez = 1000;                        

%% Geometria
% N�s [NumNo]
Nos = [1 ; 2 ; 3 ; 4];

%% Elementos
% Spring [NumElem, No'i', No'j']
Spring = [1 1 2; 2 4 2; 3 2 3];

%% Restri�oes de apoio
% Restri��es de apoio [NumNo, Rx, Ry]
Restr = [1 1 1; 3 1 1; 4 1 1];

%% Carregamento
% Carregametno [No, Px, Py]
P = [2 7.07*sqrt(2)/2 7.07*sqrt(2)/2];

%% Montagem inicial de matrizes e vetores
% Grau de liberdade
GL = size(Nos, 1) * 2;
% N�mero de elementos 
NumElem = size(Spring, 1);
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
% Coordenadas por N�
Coor = zeros(NumElem, 4);
% Seno e coseno
CS = zeros(NumElem, 2);

%% Matriz de rigidez K (Global)
for elem = 1:NumElem
    % N�s 'i' 'j'
    NoI = Spring(elem, 2);
    NoJ = Spring(elem, 3);
    % Determina��o das coordenadas nodais 
    Coor(elem, :) = [NoI * 2 - [1, 0], NoJ * 2 - [1, 0]];
    % Matriz de transforma��o T
    CS = [1 0; 0 1; sqrt(3)/2 -0.5];
    C = CS(elem, 1); S = CS(elem, 2);
    T = [C S 0 0; 0 0 C S];
    % Determina��o Ke
    Ke(:, :, elem) = Rigidez * [1  -1; -1 1];
    % Determina��o K
    K(Coor(elem, :) , Coor(elem, :)) = K(Coor(elem, :) , Coor(elem, :)) + T' *  Ke(:, :, elem) * T;
    KT(:, :, elem) = T' * (Ke(:, :, elem)) * T;
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
