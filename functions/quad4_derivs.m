function [dN, N] = quad4_derivs (r, s)
%  QUAD4_DERIVS Função de elemetno quadrilateral com 4 nós que retorna a matriz de devivadas 
% $\xi \;e\;\eta \;\;$e sua função de forma para $\xi \;e\;\eta \;\;$inseridos
% 
% dN = derivadas
% 
% N  = funão de forma
% 
% r = ponto avaliado em $\xi \;$
% 
% s = ponto avaliado em $\eta \;$
syms xi eta
n = [1.0/4.0 * (1 - xi) * (1 - eta)
     1.0/4.0 * (1 + xi) * (1 - eta)
     1.0/4.0 * (1 + xi) * (1 + eta)
     1.0/4.0 * (1 - xi) * (1 + eta)];
dN = subs([diff(n, xi) diff(n, eta)], [xi eta], [r s]);
N  = subs(n, [xi eta], [r s]);
end