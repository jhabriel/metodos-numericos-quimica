function test_edos()
% Euler explicito, Euler implicito y Heun.

% Escalar: y' = -2y, y(0)=1 -> y = exp(-2t)
f = @(t, y) -2 * y;
n = 200;
[t, Y] = euler_explicito(f, [0 1], 1, n);
eE = abs(Y(end) - exp(-2));
[~, Y] = euler_implicito(-2, [0 1], 1, n);
eI = abs(Y(end) - exp(-2));
[~, Y] = heun(f, [0 1], 1, n);
eH = abs(Y(end) - exp(-2));
assert(eE < 2e-3 && eI < 2e-3, 'Euler error O(h)');
assert(eH < 1e-4 && eH < eE / 10, 'Heun O(h^2)');

% Sistema lineal: A-> B -> C con k1=1, k2=0.5
k1 = 1; k2 = 0.5;
g = @(t, c) [-k1 * c(1); k1 * c(1) - k2 * c(2); k2 * c(2)];
[~, C] = heun(g, [0 10], [1; 0; 0], 2000);
assert(max(abs(sum(C, 2) - 1)) < 1e-10, 'conservacion de masa');
cB = @(t) k1 / (k2 - k1) * (exp(-k1 * t) - exp(-k2 * t));
assert(abs(C(end, 2) - cB(10)) < 1e-4, 'perfil de B');

[~, Ci] = euler_implicito([-k1 0 0; k1 -k2 0; 0 k2 0], [0 10], [1; 0; 0], 2000);
assert(abs(Ci(end, 2) - cB(10)) < 5e-3, 'Euler implicito, sistema lineal');

% Euler implicito recibe la matriz A, no un handle
ok = false;
try, euler_implicito(@(t, y) -y, [0 1], 1, 10); catch, ok = true; end
assert(ok, 'euler_implicito debe rechazar un handle');
ok = false;
try, euler_implicito([1 2; 3 4], [0 1], 1, 10); catch, ok = true; end
assert(ok, 'euler_implicito debe validar dimensiones');

% Rigido: y' = -1000 y, y(0) = 1; con h = 0.01 el factor por paso es
% -9 (explicito, inestable) y 1/11 (implicito, estable)
[~, Ye] = euler_explicito(@(t, y) -1000 * y, [0 1], 1, 100);
[~, Yi] = euler_implicito(-1000, [0 1], 1, 100);
assert(abs(Ye(end)) > 1e6, 'explicito inestable en problema rigido');
assert(abs(Yi(end)) < 1e-6 && Yi(end) >= 0, 'implicito estable en problema rigido');
end
