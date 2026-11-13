function [L, U, P] = lu_factorizacion(A)
% LU_FACTORIZACION  Factorizacion P*A = L*U con pivoteo parcial por filas.
%
%   [L, U, P] = lu_factorizacion(A)
%
%   L : triangular inferior con diagonal unitaria
%   U : triangular superior
%   P : matriz de permutacion

n = size(A, 1);
if size(A, 2) ~= n
    error('lu_factorizacion:dimension', 'A debe ser cuadrada.');
end

U = A;
L = eye(n);
P = eye(n);

for k = 1:n-1
    [pmax, idx] = max(abs(U(k:n, k)));
    if pmax == 0
        error('lu_factorizacion:singular', 'Matriz singular (columna %d nula).', k);
    end
    p = idx + k - 1;
    if p ~= k
        U([k p], :) = U([p k], :);
        P([k p], :) = P([p k], :);
        if k > 1
            L([k p], 1:k-1) = L([p k], 1:k-1);
        end
    end
    for i = k+1:n
        L(i, k) = U(i, k) / U(k, k);
        U(i, k:n) = U(i, k:n) - L(i, k) * U(k, k:n);
        U(i, k) = 0;
    end
end
if U(n, n) == 0
    error('lu_factorizacion:singular', 'Matriz singular.');
end
end
