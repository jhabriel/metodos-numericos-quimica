function [x, iter, historia] = biseccion(f, a, b, tol, maxit)
% BISECCION  Raiz de f(x)=0 en [a,b] por el metodo de biseccion.
%
%   [x, iter, historia] = biseccion(f, a, b, tol, maxit)
%
%   f     : handle de la funcion, f(x) escalar
%   a, b  : extremos del intervalo, con f(a)*f(b) < 0
%   tol   : tolerancia sobre el ancho del intervalo (defecto 1e-8)
%   maxit : maximo de iteraciones (defecto 100)
%
%   x     : aproximacion de la raiz
%   iter  : iteraciones realizadas
%   historia : matriz [k, a, b, c, f(c)] con la historia

if nargin < 4, tol = 1e-8; end
if nargin < 5, maxit = 100; end

fa = f(a); fb = f(b);
if fa * fb > 0
    error('biseccion:intervalo', 'f(a) y f(b) deben tener signos opuestos.');
end

historia = zeros(0, 5);
x = a;
for iter = 1:maxit
    c = (a + b) / 2;
    fc = f(c);
    historia(end+1, :) = [iter, a, b, c, fc]; %#ok<AGROW>
    if fc == 0 || (b - a) / 2 < tol
        x = c;
        return
    end
    if fa * fc < 0
        b = c; fb = fc; %#ok<NASGU>
    else
        a = c; fa = fc;
    end
    x = c;
end
warning('biseccion:maxit', 'Se alcanzo el maximo de iteraciones.');
end
