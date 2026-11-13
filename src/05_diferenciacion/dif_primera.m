function d = dif_primera(f, x, h, tipo)
% DIF_PRIMERA  Primera derivada numerica por diferencias finitas.
%
%   d = dif_primera(f, x, h, tipo)
%
%   f    : handle de la funcion (vectorizado)
%   x    : punto(s) de evaluacion
%   h    : paso
%   tipo : 'adelante' (O(h))      f'(x) ~ (f(x+h)-f(x))/h
%          'atras'    (O(h))      f'(x) ~ (f(x)-f(x-h))/h
%          'central'  (O(h^2))    f'(x) ~ (f(x+h)-f(x-h))/(2h)   [defecto]

if nargin < 4, tipo = 'central'; end

switch lower(tipo)
    case 'adelante'
        d = (f(x + h) - f(x)) / h;
    case 'atras'
        d = (f(x) - f(x - h)) / h;
    case 'central'
        d = (f(x + h) - f(x - h)) / (2 * h);
    otherwise
        error('dif_primera:tipo', 'tipo debe ser ''adelante'', ''atras'' o ''central''.');
end
end
