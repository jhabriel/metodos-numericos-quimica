function [x, iter, historia] = newton(f, df, x0, tol, maxit)
% NEWTON  Raiz de f(x)=0 por el metodo de Newton-Raphson (1 variable).
%
%   [x, iter, historia] = newton(f, df, x0, tol, maxit)
%
%   f, df : handles de f y su derivada
%   x0    : estimacion inicial
%   tol   : tolerancia sobre |x_{k+1}-x_k| (defecto 1e-10)
%   maxit : maximo de iteraciones (defecto 50)
%
%   historia : matriz [k, x_k, f(x_k)]

if nargin < 4, tol = 1e-10; end
if nargin < 5, maxit = 50; end

x = x0;
historia = zeros(0, 3);
for iter = 1:maxit
    fx = f(x);
    dfx = df(x);
    if dfx == 0
        error('newton:derivadaCero', 'Derivada nula en x = %g.', x);
    end
    xn = x - fx / dfx;
    historia(end+1, :) = [iter, xn, f(xn)]; %#ok<AGROW>
    if abs(xn - x) < tol
        x = xn;
        return
    end
    x = xn;
end
warning('newton:maxit', 'Se alcanzo el maximo de iteraciones.');
end
