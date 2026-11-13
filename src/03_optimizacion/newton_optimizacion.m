function [x, iter, hist, tipo] = newton_optimizacion(df, d2f, x0, tol, maxit)
% NEWTON_OPTIMIZACION  Punto estacionario de f (1 variable) por Newton.
%
%   [x, iter, hist, tipo] = newton_optimizacion(df, d2f, x0, tol, maxit)
%
%   Aplica Newton a f'(x) = 0:   x_{k+1} = x_k - f'(x_k) / f''(x_k)
%
%   df, d2f : handles de la primera y segunda derivada de f
%   x0      : estimacion inicial
%   tol     : tolerancia sobre |x_{k+1}-x_k| (defecto 1e-10)
%   maxit   : maximo de iteraciones (defecto 50)
%
%   tipo    : 'minimo' (f''>0), 'maximo' (f''<0) o 'indeterminado'
%   hist    : matriz [k, x_k, f'(x_k), f''(x_k)]
%
%   Converge rapido si x0 esta cerca del optimo, pero puede ir a otro
%   punto estacionario o divergir si esta lejos (a diferencia de
%   seccion_aurea, que no necesita derivadas).

if nargin < 4, tol = 1e-10; end
if nargin < 5, maxit = 50; end

x = x0;
hist = zeros(0, 4);
convergio = false;
for iter = 1:maxit
    g = df(x);
    H = d2f(x);
    if H == 0
        error('newton_optimizacion:segundaDerivadaCero', 'f''''(x) = 0 en x = %g.', x);
    end
    xn = x - g / H;
    hist(end+1, :) = [iter, xn, df(xn), d2f(xn)]; %#ok<AGROW>
    done = abs(xn - x) < tol;
    x = xn;
    if done
        convergio = true;
        break
    end
end
if ~convergio
    warning('newton_optimizacion:maxit', 'Se alcanzo el maximo de iteraciones.');
end

H = d2f(x);
if H > 0
    tipo = 'minimo';
elseif H < 0
    tipo = 'maximo';
else
    tipo = 'indeterminado';
end
end
