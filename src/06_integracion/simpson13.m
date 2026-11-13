function I = simpson13(f, a, b, n)
% SIMPSON13  Integral de f en [a,b] con Simpson 1/3 compuesto.
%
%   I = simpson13(f, a, b, n)
%
%   n : numero de subintervalos, debe ser PAR.

if n < 2 || mod(n, 2) ~= 0
    error('simpson13:n', 'n debe ser un entero par >= 2.');
end
h = (b - a) / n;
x = a + (0:n) * h;
y = f(x);
I = h / 3 * (y(1) + 4 * sum(y(2:2:n)) + 2 * sum(y(3:2:n-1)) + y(n+1));
end
