function x = gauss_pivoteo(A, b)
% GAUSS_PIVOTEO  Resuelve A*x = b por Gauss con pivoteo parcial por filas.
%
%   x = gauss_pivoteo(A, b)
%
%   En cada columna k se intercambia la fila k con la fila de mayor
%   valor absoluto en A(k:n, k).

n = size(A, 1);
b = b(:);
if size(A, 2) ~= n || length(b) ~= n
    error('gauss_pivoteo:dimension', 'A debe ser cuadrada y b de longitud compatible.');
end

for k = 1:n-1
    [pmax, idx] = max(abs(A(k:n, k)));
    if pmax == 0
        error('gauss_pivoteo:singular', 'Matriz singular (columna %d nula).', k);
    end
    p = idx + k - 1;
    if p ~= k
        A([k p], :) = A([p k], :);
        b([k p]) = b([p k]);
    end
    for i = k+1:n
        m = A(i, k) / A(k, k);
        A(i, k:n) = A(i, k:n) - m * A(k, k:n);
        b(i) = b(i) - m * b(k);
    end
end

x = sustitucion_atras(triu(A), b);
end
