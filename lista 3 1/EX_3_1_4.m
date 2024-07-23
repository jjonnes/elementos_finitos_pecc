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
% Nós [NumNo]
Nos = [1 ; 2 ; 3 ; 4];

%% Elementos
% Spring [NumElem, No'i', No'j']
Spring = [1 1 2; 2 4 2; 3 2 3];

%% Restriçoes de apoio
% Restrições de apoio [NumNo, Rx, Ry]
Restr = [1 1 1; 3 1 1; 4 1 1];

%% Carregamento
% Carregametno [No, Px, Py]
P = [2 7.07*sqrt(2)/2 7.07*sqrt(2)/2];

%% Montagem inicial de matrizes e vetores
% Grau de liberdade
GL = size(Nos, 1) * 2;
% Número de elementos 
NumElem = size(Spring, 1);
% Matriz K (elemento)
Ke = zeros(2, 2, NumElem);
% Matriz K (Global)
K = zeros(GL, GL);
% Vetor U
U = zeros(GL, 1);
% Vetor F
F = zeros(GL, 1);
% Força de reação
Fr = zeros(GL, 1);
% Coordenadas por Nó
Coor = zeros(NumElem, 4);
% Seno e coseno
CS = zeros(NumElem, 2);

%% Matriz de rigidez K (Global)
for elem = 1:NumElem
    % Nós 'i' 'j'
    NoI = Spring(elem, 2);
    NoJ = Spring(elem, 3);
    % Determinação das coordenadas nodais 
    Coor(elem, :) = [NoI * 2 - [1, 0], NoJ * 2 - [1, 0]];
    % Matriz de transformação T
    CS = [1 0; 0 1; sqrt(3)/2 -0.5];
    C = CS(elem, 1); S = CS(elem, 2);
    T = [C S 0 0; 0 0 C S];
    % Determinação Ke
    Ke(:, :, elem) = Rigidez * [1  -1; -1 1];
    % Determinação K
    K(Coor(elem, :) , Coor(elem, :)) = K(Coor(elem, :) , Coor(elem, :)) + T' *  Ke(:, :, elem) * T;
    KT(:, :, elem) = T' * (Ke(:, :, elem)) * T;
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
U(GLSR) = K(GLSR, GLSR) \ F(GLSR);

%% Forças de reação
% Força de reação
Fr(GLR) = K(GLR, :) * U;
