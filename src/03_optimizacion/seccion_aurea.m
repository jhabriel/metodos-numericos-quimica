function [xopt, fopt, iter] = seccion_aurea(f, a, b, tol, tipo, maxit)
% SECCION_AUREA  Optimo de f en [a,b] (unimodal) por seccion dorada.
%
%   [xopt, fopt, iter] = seccion_aurea(f, a, b, tol, tipo, maxit)
%
%   f     : handle de la funcion, unimodal en [a,b]
%   tol   : tolerancia sobre el ancho del intervalo (defecto 1e-6)
%   tipo  : 'min' (defecto) o 'max'
%   maxit : maximo de iteraciones (defecto 200)

if nargin < 4, tol = 1e-6; end
if nargin < 5, tipo = 'min'; end
if nargin < 6, maxit = 200; end

signo = 1;
if strcmp(tipo, 'max'), signo = -1; end
g = @(x) signo * f(x);

R = (sqrt(5) - 1) / 2;          % 0.6180...
x1 = b - R * (b - a);
x2 = a + R * (b - a);
g1 = g(x1); g2 = g(x2);

for iter = 1:maxit
    if (b - a) < tol, break; end
    if g1 < g2
        b = x2;
        x2 = x1; g2 = g1;
        x1 = b - R * (b - a);
        g1 = g(x1);
    else
        a = x1;
        x1 = x2; g1 = g2;
        x2 = a + R * (b - a);
        g2 = g(x2);
    end
end

xopt = (a + b) / 2;
fopt = f(xopt);
end
