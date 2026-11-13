function test_interpolacion_regresion()
% Lagrange y regresion lineal generalizada.

% Lagrange reproduce exactamente un polinomio de grado <= n-1
p = @(x) 3 * x.^3 - 2 * x + 5;
xn = [-1 0 1 2];
xq = linspace(-1, 2, 7);
assert(max(abs(lagrange(xn, p(xn), xq) - p(xq))) < 1e-12, 'lagrange polinomio');
assert(abs(lagrange(xn, p(xn), xn(2)) - p(xn(2))) < 1e-14, 'lagrange en nodo');

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

% Arrhenius: ln k = ln A - Ea/(R T)
R = 8.314; A0 = 1e10; Ea = 75e3;
T = [300 310 320 330 340]';
k = A0 * exp(-Ea ./ (R * T));
a = regresion_lineal_general(1 ./ T, log(k), {@(x) ones(size(x)), @(x) x});
assert(abs(-a(2) * R - Ea) / Ea < 1e-8, 'Arrhenius Ea');
end
