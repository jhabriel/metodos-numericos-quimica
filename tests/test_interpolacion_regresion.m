function test_interpolacion_regresion()
% Lagrange y regresion lineal generalizada.

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

% Recta exacta
x = (0:5)'; y = 2 * x + 1;
[a, r2] = regresion_lineal_general(x, y, {@(x) ones(size(x)), @(x) x});
assert(norm(a - [1; 2]) < 1e-10 && abs(r2 - 1) < 1e-12, 'recta exacta');

% Parabola con ruido: comparar con polyfit
rand('seed', 2); %#ok<RAND>
x = linspace(0, 4, 20)'; y = 1 - 2 * x + 0.5 * x.^2 + 0.01 * rand(20, 1);
a = regresion_lineal_general(x, y, {@(x) ones(size(x)), @(x) x, @(x) x.^2});
pf = polyfit(x, y, 2);
assert(norm(a - fliplr(pf)') < 1e-8, 'parabola vs polyfit');

% Regresion multiple: plano exacto y = 1 + 2*x1 - 3*x2
[X1, X2] = meshgrid(0:2, 0:3);
X = [X1(:), X2(:)]; y = 1 + 2 * X(:, 1) - 3 * X(:, 2);
bm = {@(X) ones(size(X, 1), 1), @(X) X(:, 1), @(X) X(:, 2)};
[a, r2] = regresion_lineal_general(X, y, bm);
assert(norm(a - [1; 2; -3]) < 1e-10 && abs(r2 - 1) < 1e-12, 'regresion multiple');

% Arrhenius: ln k = ln A - Ea/(R T)
R = 8.314; A0 = 1e10; Ea = 75e3;
T = [300 310 320 330 340]';
k = A0 * exp(-Ea ./ (R * T));
a = regresion_lineal_general(1 ./ T, log(k), {@(x) ones(size(x)), @(x) x});
assert(abs(-a(2) * R - Ea) / Ea < 1e-8, 'Arrhenius Ea');
end
