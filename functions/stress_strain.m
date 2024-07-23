function [sig, eps] = stress_strain(nodes,...
    element,displacement,D,q,coordinateelem)

nodeselement= size(element, 2);
ndof        = size(nodes, 2);
numelement  = size(element, 1);
sig         = zeros(numelement, q, 3);
for i = 1: numelement
    C = [nodes(element(i,:), :)];
    U = displacement(coordinateelem(i,:));
    for j = 1: q
         B = compute_B(C, q, nodeselement, ndof);
         sig(i,j,:) = (D*B*U);
         eps = B*U;
    end
end
end