function yq = lagrange(x, y, xq)
% LAGRANGE  Evalua en xq el polinomio interpolante de Lagrange.
%
%   yq = lagrange(x, y, xq)
%
%   x, y : nodos (n puntos, x distintos) -> polinomio de grado n-1
%   xq   : puntos de evaluacion (escalar o arreglo)

x = x(:); y = y(:);
n = length(x);
if length(y) ~= n
    error('lagrange:dimension', 'x e y deben tener la misma longitud.');
end
if numel(unique(x)) ~= n
    error('lagrange:nodosRepetidos', 'Los nodos x deben ser distintos.');
end

yq = zeros(size(xq));
for i = 1:n
    Li = ones(size(xq));
    for j = [1:i-1, i+1:n]
        Li = Li .* (xq - x(j)) / (x(i) - x(j));
    end
    yq = yq + y(i) * Li;
end
end
