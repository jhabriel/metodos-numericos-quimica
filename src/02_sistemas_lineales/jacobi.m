function [x, iter, hist] = jacobi(A, b, x0, tol, maxit)
% JACOBI  Resuelve A*x = b por el metodo iterativo de Jacobi.
%
%   [x, iter, hist] = jacobi(A, b, x0, tol, maxit)
%
%   A     : matriz cuadrada con diagonal no nula (convergencia garantizada
%           si es diagonalmente dominante por filas)
%   b     : lado derecho
%   x0    : estimacion inicial (defecto: ceros)
%   tol   : tolerancia sobre norm(x_{k+1}-x_k, inf) (defecto 1e-8)
%   maxit : maximo de iteraciones (defecto 500)
%
%   x_i^{(k+1)} = ( b_i - sum_{j~=i} a_ij x_j^{(k)} ) / a_ii
%
%   hist  : matriz [k, norm(dx, inf)]

n = size(A, 1);
b = b(:);
if size(A, 2) ~= n || length(b) ~= n
    error('jacobi:dimension', 'A debe ser cuadrada y b de longitud compatible.');
end
if nargin < 3 || isempty(x0), x0 = zeros(n, 1); end
if nargin < 4, tol = 1e-8; end
if nargin < 5, maxit = 500; end
if any(diag(A) == 0)
    error('jacobi:diagonalNula', 'A tiene elementos nulos en la diagonal.');
end

x = x0(:);
hist = zeros(0, 2);
for iter = 1:maxit
    xn = zeros(n, 1);
    for i = 1:n
        s = A(i, :) * x - A(i, i) * x(i);     % suma sobre j ~= i con x viejo
        xn(i) = (b(i) - s) / A(i, i);
    end
    err = norm(xn - x, inf);
    hist(end+1, :) = [iter, err]; %#ok<AGROW>
    x = xn;
    if err < tol
        return
    end
end
warning('jacobi:maxit', 'Se alcanzo el maximo de iteraciones sin converger.');
end
