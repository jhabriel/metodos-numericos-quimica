function I = trapecio(f, a, b, n)
% TRAPECIO  Integral de f en [a,b] con la regla del trapecio compuesta.
%
%   I = trapecio(f, a, b, n)
%
%   f : handle (vectorizado)
%   n : numero de subintervalos (n >= 1)

if n < 1 || n ~= round(n)
    error('trapecio:n', 'n debe ser un entero positivo.');
end
h = (b - a) / n;
x = a + (0:n) * h;
y = f(x);
I = h / 2 * (y(1) + 2 * sum(y(2:end-1)) + y(end));
end
