function test_interpolacion_regresion()
% Lagrange, Vandermonde y regresion por ecuaciones normales.

% Lagrange: coeficientes de un polinomio de grado <= n-1 (potencias descendentes)
p = @(x) 3 * x.^3 - 2 * x + 5;
xn = [-1 0 1 2];
xq = linspace(-1, 2, 7);
[c, L] = lagrange(xn, p(xn));
assert(norm(c - [3 0 -2 5]) < 1e-10, 'lagrange coeficientes');
assert(max(abs(polyval(c, xq) - p(xq))) < 1e-10, 'lagrange evaluado');
assert(norm(c - p(xn) * L) < 1e-12, 'p = y''*L');
for i = 1:4          % la base vale 1 en su nodo y 0 en los demas
    for j = 1:4
        assert(abs(polyval(L(i, :), xn(j)) - (i == j)) < 1e-12, 'base de Lagrange');
    end
end
assert(norm(lagrange([0 1 2], [1 3 7]) - [1 1 1]) < 1e-12, 'lagrange x^2+x+1');
ok = false;
try, lagrange([1 1 2], [1 2 3]); catch, ok = true; end
assert(ok, 'lagrange debe rechazar nodos repetidos');

% Vandermonde armada a mano: mismo polinomio que Lagrange
xv = [0 1 2 4]';
yv = [1 3 7 21]';
V = [xv.^3, xv.^2, xv, ones(size(xv))];
pv = resolver_lu(V, yv)';
assert(norm(pv - lagrange(xv, yv)) < 1e-10, 'vandermonde = lagrange');
assert(norm(V * pv' - yv) < 1e-10, 'V*p = y');

% Regresion por ecuaciones normales (Z'Z) a = Z'y, resueltas con resolver_lu
% Recta exacta
x = (0:5)'; y = 2 * x + 1;
Z = [ones(size(x)), x];
a = resolver_lu(Z' * Z, Z' * y);
assert(norm(a - [1; 2]) < 1e-10, 'recta exacta');

% Parabola con ruido: comparar con polyfit y con A\b
rand('seed', 2); %#ok<RAND>
x = linspace(0, 4, 20)'; y = 1 - 2 * x + 0.5 * x.^2 + 0.01 * rand(20, 1);
Z = [ones(size(x)), x, x.^2];
a = resolver_lu(Z' * Z, Z' * y);
pf = polyfit(x, y, 2);
assert(norm(a - fliplr(pf)') < 1e-8, 'parabola vs polyfit');
assert(norm(a - (Z' * Z) \ (Z' * y)) < 1e-10, 'parabola vs A\b');

% Regresion multiple: plano exacto y = 1 + 2*x1 - 3*x2
[X1, X2] = meshgrid(0:2, 0:3);
x1 = X1(:); x2 = X2(:); y = 1 + 2 * x1 - 3 * x2;
Z = [ones(size(x1)), x1, x2];
a = resolver_lu(Z' * Z, Z' * y);
assert(norm(a - [1; 2; -3]) < 1e-10, 'regresion multiple');

% Arrhenius: ln k = ln A - Ea/(R T)
R = 8.314; A0 = 1e10; Ea = 75e3;
T = [300 310 320 330 340]';
k = A0 * exp(-Ea ./ (R * T));
Z = [ones(size(T)), 1 ./ T];
a = resolver_lu(Z' * Z, Z' * log(k));
assert(abs(-a(2) * R - Ea) / Ea < 1e-8, 'Arrhenius Ea');
end
