function [x, iter, hist] = newton_sistema(F, J, x0, tol, maxit)
% NEWTON_SISTEMA  Resuelve F(x)=0 (sistema no lineal) con Newton-Raphson.
%
%   [x, iter, hist] = newton_sistema(F, J, x0, tol, maxit)
%
%   F     : handle, F(x) devuelve vector columna n x 1
%   J     : handle del jacobiano n x n; si es [] se aproxima por
%           diferencias finitas centradas
%   x0    : estimacion inicial (vector columna)
%   tol   : tolerancia sobre norm(dx) (defecto 1e-10)
%   maxit : maximo de iteraciones (defecto 50)
%
%   En cada paso se resuelve J*dx = -F con factorizacion LU propia.
%   hist  : matriz [k, norm(dx), norm(F)]

if nargin < 4, tol = 1e-10; end
if nargin < 5, maxit = 50; end

x = x0(:);
hist = zeros(0, 3);
for iter = 1:maxit
    Fx = F(x);
    Fx = Fx(:);
    if isempty(J)
        Jx = jacobiano_numerico(F, x);
    else
        Jx = J(x);
    end
    dx = resolver_lu(Jx, -Fx);
    x = x + dx;
    hist(end+1, :) = [iter, norm(dx), norm(F(x))]; %#ok<AGROW>
    if norm(dx) < tol
        return
    end
end
warning('newton_sistema:maxit', 'Se alcanzo el maximo de iteraciones.');
end

function Jx = jacobiano_numerico(F, x)
n = numel(x);
Jx = zeros(n);
for j = 1:n
    h = 1e-6 * max(1, abs(x(j)));
    e = zeros(n, 1); e(j) = h;
    Fp = F(x + e);
    Fm = F(x - e);
    Jx(:, j) = (Fp(:) - Fm(:)) / (2 * h);
end
end
