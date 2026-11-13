function test_optimizacion()
% Seccion aurea.

[x, fx] = seccion_aurea(@(x) (x - 2).^2 + 1, 0, 5, 1e-8);
assert(abs(x - 2) < 1e-6 && abs(fx - 1) < 1e-10, 'minimo');

[x, fx] = seccion_aurea(@(x) 2 * sin(x) - x.^2 / 10, 0, 4, 1e-8, 'max');
assert(abs(x - 1.427552) < 1e-5 && abs(fx - 1.775726) < 1e-5, 'maximo');
end
