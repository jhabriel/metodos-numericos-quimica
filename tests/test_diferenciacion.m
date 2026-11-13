function test_diferenciacion()
% Derivadas numericas y sus ordenes de convergencia.

f = @(x) exp(x); x0 = 1; ex = exp(1);
h = 1e-3;

dadl = dif_primera(f, x0, h, 'adelante');
datr = dif_primera(f, x0, h, 'atras');
dcen = dif_primera(f, x0, h, 'central');
assert(abs(dadl - ex) < 2e-3 && abs(datr - ex) < 2e-3, 'adelante/atras O(h)');
assert(abs(dcen - ex) < 1e-5, 'central O(h^2)');
assert(abs(dcen - ex) < abs(dadl - ex), 'central mejor que adelante');

% Orden: reducir h a la mitad -> error / 2 (O(h)) o / 4 (O(h^2))
e1 = abs(dif_primera(f, x0, 0.01, 'adelante') - ex);
e2 = abs(dif_primera(f, x0, 0.005, 'adelante') - ex);
assert(abs(e1 / e2 - 2) < 0.1, 'orden 1');
e1 = abs(dif_primera(f, x0, 0.01, 'central') - ex);
e2 = abs(dif_primera(f, x0, 0.005, 'central') - ex);
assert(abs(e1 / e2 - 4) < 0.1, 'orden 2');

d2 = dif_segunda_central(f, x0, 1e-3);
assert(abs(d2 - ex) < 1e-5, 'segunda derivada');

% Exacta para cuadraticas (segunda) y lineales (central)
assert(abs(dif_segunda_central(@(x) 3 * x.^2, 5, 0.1) - 6) < 1e-10, 'cuadratica');
end
