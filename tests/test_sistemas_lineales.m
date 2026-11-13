function test_sistemas_lineales()
% Sustituciones, Gauss, LU.

A = [4 -2 1; -2 4 -2; 1 -2 4];
b = [11; -16; 17];
xref = A \ b;

assert(norm(gauss_naive(A, b) - xref) < 1e-12, 'gauss_naive');
assert(norm(gauss_pivoteo(A, b) - xref) < 1e-12, 'gauss_pivoteo');
assert(norm(resolver_lu(A, b) - xref) < 1e-12, 'resolver_lu');

[L, U, P] = lu_factorizacion(A);
assert(norm(P * A - L * U) < 1e-12, 'P*A = L*U');
assert(isequal(L, tril(L)) && isequal(U, triu(U)), 'forma de L y U');

% Matriz que exige pivoteo (pivote nulo en (1,1))
B = [0 2 1; 1 1 1; 2 1 3];
c = [3; 3; 7];
assert(norm(gauss_pivoteo(B, c) - B \ c) < 1e-12, 'pivoteo con a11 = 0');
assert(norm(resolver_lu(B, c) - B \ c) < 1e-12, 'LU con a11 = 0');
ok = false;
try, gauss_naive(B, c); catch, ok = true; end
assert(ok, 'gauss_naive debe fallar con pivote nulo');

% Sustituciones
L2 = [2 0 0; 1 3 0; 4 5 6];
assert(norm(sustitucion_adelante(L2, [2; 7; 32]) - L2 \ [2; 7; 32]) < 1e-12, 'adelante');
U2 = L2';
assert(norm(sustitucion_atras(U2, [2; 7; 32]) - U2 \ [2; 7; 32]) < 1e-12, 'atras');

% Metodos iterativos (diagonalmente dominante)
D = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8];
d = [6; 25; -11; 15];
[xj, itj] = jacobi(D, d, [], 1e-12);
[xg, itg] = gauss_seidel(D, d, [], 1e-12);
assert(norm(xj - [1; 2; -1; 1]) < 1e-10, 'jacobi');
assert(norm(xg - [1; 2; -1; 1]) < 1e-10, 'gauss_seidel');
assert(itg < itj, 'Gauss-Seidel converge en menos iteraciones que Jacobi aqui');
ok = false;
try, jacobi([0 1; 1 0], [1; 1]); catch, ok = true; end
assert(ok, 'jacobi debe rechazar diagonal nula');

% Aleatoria
rand('seed', 1); %#ok<RAND>
M = rand(6) + 6 * eye(6); r = rand(6, 1);
assert(norm(resolver_lu(M, r) - M \ r) < 1e-10, 'LU aleatoria');
end
