function test_no_lineales()
% Biseccion, Newton y Newton para sistemas.

f  = @(x) x.^2 - 2;
df = @(x) 2 * x;

x = biseccion(f, 0, 2, 1e-10);
assert(abs(x - sqrt(2)) < 1e-9, 'biseccion');

x = newton(f, df, 1, 1e-12);
assert(abs(x - sqrt(2)) < 1e-12, 'newton');

% Sistema: x^2 + y^2 = 4,  x*y = 1
F = @(v) [v(1)^2 + v(2)^2 - 4; v(1) * v(2) - 1];
J = @(v) [2*v(1), 2*v(2); v(2), v(1)];
v = newton_sistema(F, J, [2; 0.5]);
assert(norm(F(v)) < 1e-10, 'newton_sistema con jacobiano');

v2 = newton_sistema(F, [], [2; 0.5]);   % jacobiano numerico
assert(norm(F(v2)) < 1e-8, 'newton_sistema jacobiano numerico');

% Intervalo sin cambio de signo debe fallar
ok = false;
try, biseccion(f, 2, 3); catch, ok = true; end
assert(ok, 'biseccion debe rechazar intervalo sin cambio de signo');
end
