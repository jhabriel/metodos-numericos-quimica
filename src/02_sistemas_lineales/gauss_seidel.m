function [x, iter, historia] = gauss_seidel(A, b, x0, tol, maxit)
% GAUSS_SEIDEL  Resuelve A*x = b por el metodo iterativo de Gauss-Seidel.
%
%   [x, iter, historia] = gauss_seidel(A, b, x0, tol, maxit)
%
%   Mismos argumentos que jacobi. La diferencia es que cada componente
%   usa de inmediato los valores ya actualizados en la misma iteracion:
%
%   x_i^{(k+1)} = ( b_i - sum_{j<i} a_ij x_j^{(k+1)} - sum_{j>i} a_ij x_j^{(k)} ) / a_ii
%
%   historia : matriz [k, norm(dx, inf)]

n = size(A, 1);
b = b(:);
if size(A, 2) ~= n || length(b) ~= n
    error('gauss_seidel:dimension', 'A debe ser cuadrada y b de longitud compatible.');
end
if nargin < 3 || isempty(x0), x0 = zeros(n, 1); end
if nargin < 4, tol = 1e-8; end
if nargin < 5, maxit = 500; end
if any(diag(A) == 0)
    error('gauss_seidel:diagonalNula', 'A tiene elementos nulos en la diagonal.');
end

x = x0(:);
historia = zeros(0, 2);
for iter = 1:maxit
    xold = x;
    for i = 1:n
        s = A(i, :) * x - A(i, i) * x(i);     % x ya contiene los valores nuevos j < i
        x(i) = (b(i) - s) / A(i, i);
    end
    err = norm(x - xold, inf);
    historia(end+1, :) = [iter, err]; %#ok<AGROW>
    if err < tol
        return
    end
end
warning('gauss_seidel:maxit', 'Se alcanzo el maximo de iteraciones sin converger.');
end
