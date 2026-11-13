function I = simpson38(f, a, b, n)
% SIMPSON38  Integral de f en [a,b] con Simpson 3/8 compuesto.
%
%   I = simpson38(f, a, b, n)
%
%   n : numero de subintervalos, debe ser MULTIPLO DE 3.

if n < 3 || mod(n, 3) ~= 0
    error('simpson38:n', 'n debe ser un multiplo de 3 (>= 3).');
end
h = (b - a) / n;
x = a + (0:n) * h;
y = f(x);
idx = 2:n;                         % nodos interiores (indices base 1)
mult3 = mod(idx - 1, 3) == 0;      % nodos x_3, x_6, ... pesan 2
I = 3 * h / 8 * (y(1) + y(n+1) + 3 * sum(y(idx(~mult3))) + 2 * sum(y(idx(mult3))));
end
