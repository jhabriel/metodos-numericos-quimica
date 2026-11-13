function I = boole(f, a, b, n)
% BOOLE  Integral de f en [a,b] con la regla de Boole compuesta.
%
%   I = boole(f, a, b, n)
%
%   n : numero de subintervalos, debe ser MULTIPLO DE 4.
%   Cada panel de 4 subintervalos usa:
%     2h/45 * (7 f0 + 32 f1 + 12 f2 + 32 f3 + 7 f4)

if n < 4 || mod(n, 4) ~= 0
    error('boole:n', 'n debe ser un multiplo de 4 (>= 4).');
end
h = (b - a) / n;
x = a + (0:n) * h;
y = f(x);
I = 0;
for k = 1:4:n
    I = I + 2 * h / 45 * (7 * y(k) + 32 * y(k+1) + 12 * y(k+2) + 32 * y(k+3) + 7 * y(k+4));
end
end
