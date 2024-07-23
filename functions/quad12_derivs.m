function [dN, N] = quad12_derivs (r, s)
% Função de elemetno quadrilateral com 12 nós que retorna a matriz de devivadas e sua função de forma para inseridos
% dN = derivadas
% [dN, N] = quad8_derivs (xi, eta)
% N  = funão de forma
% r = ponto avaliado em xi
% s = ponto avaliado em eta
% Para resultado simbílico, chamar syms xi eta antes da função
syms xi eta
n = [
(5*(eta))/16 + (5*(xi))/16 - (5*(eta)*(xi))/16 + (9*(eta)^2)/32 - (9*(eta)^3)/32 + (9*(xi)^2)/32 - (9*(xi)^3)/32 - (9*(eta)*(xi)^2)/32 - (9*(eta)^2*(xi))/32 + (9*(eta)*(xi)^3)/32 + (9*(eta)^3*(xi))/32 - 5/16
(5*(eta))/16 - (5*(xi))/16 + (5*(eta)*(xi))/16 + (9*(eta)^2)/32 - (9*(eta)^3)/32 + (9*(xi)^2)/32 + (9*(xi)^3)/32 - (9*(eta)*(xi)^2)/32 + (9*(eta)^2*(xi))/32 - (9*(eta)*(xi)^3)/32 - (9*(eta)^3*(xi))/32 - 5/16
(9*(eta)^2)/32 - (5*(xi))/16 - (5*(eta)*(xi))/16 - (5*(eta))/16 + (9*(eta)^3)/32 + (9*(xi)^2)/32 + (9*(xi)^3)/32 + (9*(eta)*(xi)^2)/32 + (9*(eta)^2*(xi))/32 + (9*(eta)*(xi)^3)/32 + (9*(eta)^3*(xi))/32 - 5/16
(5*(xi))/16 - (5*(eta))/16 + (5*(eta)*(xi))/16 + (9*(eta)^2)/32 + (9*(eta)^3)/32 + (9*(xi)^2)/32 - (9*(xi)^3)/32 + (9*(eta)*(xi)^2)/32 - (9*(eta)^2*(xi))/32 - (9*(eta)*(xi)^3)/32 - (9*(eta)^3*(xi))/32 - 5/16
                                                                                                  (27*(eta)*(xi))/32 - (27*(xi))/32 - (9*(eta))/32 - (9*(xi)^2)/32 + (27*(xi)^3)/32 + (9*(eta)*(xi)^2)/32 - (27*(eta)*(xi)^3)/32 + 9/32
                                                                                                  (27*(xi))/32 - (9*(eta))/32 - (27*(eta)*(xi))/32 - (9*(xi)^2)/32 - (27*(xi)^3)/32 + (9*(eta)*(xi)^2)/32 + (27*(eta)*(xi)^3)/32 + 9/32
                                                                                                (9*(xi))/32 - (27*(eta))/32 - (27*(eta)*(xi))/32 - (9*(eta)^2)/32 + (27*(eta)^3)/32 - (9*(eta)^2*(xi))/32 + (27*(eta)^3*(xi))/32 + 9/32
                                                                                                (27*(eta))/32 + (9*(xi))/32 + (27*(eta)*(xi))/32 - (9*(eta)^2)/32 - (27*(eta)^3)/32 - (9*(eta)^2*(xi))/32 - (27*(eta)^3*(xi))/32 + 9/32
                                                                                                  (9*(eta))/32 + (27*(xi))/32 + (27*(eta)*(xi))/32 - (9*(xi)^2)/32 - (27*(xi)^3)/32 - (9*(eta)*(xi)^2)/32 - (27*(eta)*(xi)^3)/32 + 9/32
                                                                                                  (9*(eta))/32 - (27*(xi))/32 - (27*(eta)*(xi))/32 - (9*(xi)^2)/32 + (27*(xi)^3)/32 - (9*(eta)*(xi)^2)/32 + (27*(eta)*(xi)^3)/32 + 9/32
                                                                                                (27*(eta))/32 - (9*(xi))/32 - (27*(eta)*(xi))/32 - (9*(eta)^2)/32 - (27*(eta)^3)/32 + (9*(eta)^2*(xi))/32 + (27*(eta)^3*(xi))/32 + 9/32
                                                                                                (27*(eta)*(xi))/32 - (9*(xi))/32 - (27*(eta))/32 - (9*(eta)^2)/32 + (27*(eta)^3)/32 + (9*(eta)^2*(xi))/32 - (27*(eta)^3*(xi))/32 + 9/32];
dN = subs([diff(n, xi) diff(n, eta)], [xi eta], [r s]);
N  = subs(n, [xi eta], [r s]);
end