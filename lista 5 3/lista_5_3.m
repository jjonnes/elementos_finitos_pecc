% Exercises 5.3
% 1. Using the general Hooke’s law, demonstrate the relation between Young’s and bulk modulus.
% σ = (E/(1+ν)(1-2ν))(εx + εy + εz)I
% Comparando essa equação com a equação de estresse em termos do módulo volumétrico (K), temos: 
% K = (E/(3(1-2ν))) 
% Portanto, podemos concluir que a relação entre o módulo de elasticidade (E) e o módulo volumétrico (K) para um material elástico linear e isotrópico é:
% E = 3K(1-2ν)/(1+ν)
% 
% 2. Demonstrate the relation below between elastic constants G = E/2(1+ν) 
% E = 3K(1-2ν)/(1+ν)
% G = (3K(1-2ν)/(1+ν))/ 2(1+ν)
% G = 3K(1-2ν)/2(1+ν)²
% 
% 3. A cylindrical sample of concrete, with 20 cm height and 10 cm diameter, was tested in a uniaxial compression test where axial and radial deformations were registered. At a point during the elastic regime, the following values were logged: axial stress: 20 Mpa, axial deformation: 0.2 mm and radial deformation: 0.015 mm. Find the corresponding material D matrix tobe used in a plane stress linear-elastic analysis.
clear all
clc

syms E nu

d = E/((1+nu)*(1-2*nu))*[1-nu nu nu 0;
                        nu 1-nu nu 0;
                        nu nu 1-nu 0;
                        0  0  0  (1-2*nu)/2];
sig = [0.0; 20; 0.0; 0.0];
eps = [0.015/100; -0.2/200; 0.0; (0.015/100+0.2/200)];

A1 = d\sig == eps;
C = solve (A1(2));
A2 = subs(A1(1), E, C);
B = solve (A2);
D = C/(1-B^2)*[1 B 0;
                B 1 0;
                0 0 (1-B)/2];
