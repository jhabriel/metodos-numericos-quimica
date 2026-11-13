function test_optimizacion()
% Seccion aurea y Newton para optimizacion.

[x, fx] = seccion_aurea(@(x) (x - 2).^2 + 1, 0, 5, 1e-8);
assert(abs(x - 2) < 1e-6 && abs(fx - 1) < 1e-10, 'minimo');

% Maximo de f: se minimiza g = -f y se cambia el signo del valor
[x, gx] = seccion_aurea(@(x) -(2 * sin(x) - x.^2 / 10), 0, 4, 1e-8);
assert(abs(x - 1.427552) < 1e-5 && abs(-gx - 1.775726) < 1e-5, 'maximo via -f');

% Newton para optimizacion = newton aplicado a g' (con g'' como derivada)
dg  = @(x) -2 * cos(x) + x / 5;
d2g = @(x) 2 * sin(x) + 1 / 5;
xn = newton(dg, d2g, 2.5);
assert(abs(xn - 1.427552) < 1e-5 && d2g(xn) > 0, 'newton minimo de g (g'''' > 0)');

xn = newton(@(x) 2 * (x - 2), @(x) 2 + 0 * x, 10);
assert(abs(xn - 2) < 1e-10, 'newton minimo de una cuadratica');
end
