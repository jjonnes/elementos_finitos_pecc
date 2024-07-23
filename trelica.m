%% Dados de entrada
% Dados da estrutura
E = 2000;    %MN/m²                   
A = 1;       %m²                 

% Geometria
% Nós [NumNo, Cx, Cy] [m]
Nos = [1 0.8 0.6; 2 0 0.6; 3 0 0];

% Elementos
% Barra [NumElem, No'i', No'j']
Barra = [1 3 1; 2 2 1; 3 3 2];

% Restriçoes de apoio
% Restrições de apoio [NumNo, Rx, Ry]
Restr = [2 1 1; 3 1 0];

% Carregamento
% Carregametno [No, Px, Py]
P = [1 0 -10];

%% Montagem inicial de matrizes e vetores
GL = size(Nos, 1) * 2;      % Grau de liberdade
NumElem = size(Barra, 1);   % Número de elementos 
Ke = zeros(2, 2, NumElem);  % Matriz K (elemento)
K = zeros(GL, GL);          % Matriz K (Global)
U = zeros(GL, 1);           % Vetor U
F = zeros(GL, 1);           % Vetor F
Fr = zeros(GL, 1);          % Força de reação
L = zeros(NumElem, 1);      % Comprimento do elemnto
Coor = zeros(NumElem, 4);   % Coordenadas por Nó
CS = zeros(NumElem,2);      % Seno e coseno

%% Matriz de rigidez K (Global)
for elem = 1:NumElem
    NoI = Barra(elem, 2);                                   % Nós 'i' 'j'
    NoJ = Barra(elem, 3);
    L(elem) = norm(Nos(NoJ, [2, 3]) - Nos(NoI, [2, 3]));    % Determinação L
    Coor(elem, :) = [NoI * 2 - [1, 0], NoJ * 2 - [1, 0]];   % Determinação das coordenadas nodais 
    CS(elem, 1) = (Nos(NoJ, 2) - Nos(NoI, 2)) /L(elem); 
    CS(elem, 2) = (Nos(NoJ, 3) - Nos(NoI, 3)) /L(elem);
    C = CS(elem, 1); S = CS(elem, 2);
    T = [C S 0 0; 0 0 C S];                                 % Matriz de transformação T
    Ke(:, :, elem) = (E * A / L(elem)) * [1  -1; -1 1];     % Determinação Ke
    K(Coor(elem, :) , Coor(elem, :)) = K(Coor(elem, :)...
        , Coor(elem, :)) + T' *  Ke(:, :, elem) * T;        % Determinação K
    KT(:, :, elem) = T' * (Ke(:, :, elem)) * T;
end

%% Vetor de forças
for load = 1: size(P, 1)
    NoP = P(load, 1);               % Nós
    CoorP = (NoP * 2 - [1, 0])';    % Determinação das cordenadas 
    Fa = P(load, [2, 3])';          % Força aplicada
    F(CoorP) = Fa;                  % Vetor de forças      
end

%% Restrições nos apoios
NumGLR = sum(sum(Restr(:, [2, 3])));% Determinação do número de restrições
GLR = zeros(NumGLR, 1);             % Vetor deslocamentos restringidos 

%% Cordenadas com restiçoes
i = 1;
for apoio = 1:(size(Restr, 1))
    NoApoio = Restr(apoio, 1);          % Nó de apoio
    if Restr(apoio, 2) == 1
        GLR(i, 1) = NoApoio * 2 - 1;    % Restrição em x
        i = i + 1;
    end
    if Restr(apoio, 3) == 1
        GLR(i, 1) = NoApoio * 2;        % Restrição em y
        i = i + 1;
    end
end

%% Cordenadas sem restrição
GLSR = setxor((1:GL)', GLR);        % Grau de liberdade sem restrição

%% Deslocamentos sem restição
U(GLSR) = K(GLSR, GLSR) \ F(GLSR);  % Deslocamentos

%% Forças de reação
Fr(GLR) = K(GLR, :) * U;            % Força de reação

%% Resuldatos 
R = [Fr, U];
fprintf('%6s %12s\r\n','F','U');
fprintf('%6f %12f\n',U,Fr);