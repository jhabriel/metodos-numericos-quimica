function test_optimizacion()
% Seccion aurea.

[x, fx] = seccion_aurea(@(x) (x - 2).^2 + 1, 0, 5, 1e-8);
assert(abs(x - 2) < 1e-6 && abs(fx - 1) < 1e-10, 'minimo');

[x, fx] = seccion_aurea(@(x) 2 * sin(x) - x.^2 / 10, 0, 4, 1e-8, 'max');
assert(abs(x - 1.427552) < 1e-5 && abs(fx - 1.775726) < 1e-5, 'maximo');

% Newton para optimizacion: mismo problema, f' y f'' analiticas
df  = @(x) 2 * cos(x) - x / 5;
d2f = @(x) -2 * sin(x) - 1 / 5;
[xn, ~, ~, tipo] = newton_optimizacion(df, d2f, 2.5);
assert(abs(xn - 1.427552) < 1e-5 && strcmp(tipo, 'maximo'), 'newton_optimizacion maximo');

[xn, ~, ~, tipo] = newton_optimizacion(@(x) 2 * (x - 2), @(x) 2 + 0 * x, 10);
assert(abs(xn - 2) < 1e-10 && strcmp(tipo, 'minimo'), 'newton_optimizacion cuadratica');
end
