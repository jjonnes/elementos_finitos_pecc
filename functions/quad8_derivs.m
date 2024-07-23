function [dN, N] = quad8_derivs (r, s)
% Função de elemetno quadrilateral com 8 nós que retorna a matriz de devivadas e sua função de forma para inseridos
% dN = derivadas
% [dN, N] = quad8_derivs (xi, eta)
% N  = funão de forma
% r = ponto avaliado em xi
% s = ponto avaliado em eta
% Para resultado simbílico, chamar syms xi eta antes da função
syms xi eta
rp1=1.0+xi; rm1=1.0-xi;sp1=1.0+eta; sm1=1.0-eta;
n = [0.25*rm1*sm1*(rm1+sm1-3.0)
    0.25*rp1*sm1*(rp1+sm1-3.0)
    0.25*rp1*sp1*(rp1+sp1-3.0)
    0.25*rm1*sp1*(rm1+sp1-3.0)
    0.50*sm1*(1.0-xi*xi)
    0.50*rp1*(1.0-eta*eta)
    0.50*sp1*(1.0-xi*xi)
    0.50*rm1*(1.0-eta*eta)];
dN = subs([diff(n, xi) diff(n, eta)], [xi eta], [r s]);
N  = subs(n, [xi eta], [r s]);
end