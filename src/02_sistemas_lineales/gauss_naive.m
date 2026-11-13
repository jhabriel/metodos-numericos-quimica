function x = gauss_naive(A, b)
% GAUSS_NAIVE  Resuelve A*x = b por eliminacion de Gauss sin pivoteo.
%
%   x = gauss_naive(A, b)
%
%   Falla si aparece un pivote nulo; puede ser numericamente inestable
%   con pivotes pequenos (ver gauss_pivoteo).

n = size(A, 1);
b = b(:);
if size(A, 2) ~= n || length(b) ~= n
    error('gauss_naive:dimension', 'A debe ser cuadrada y b de longitud compatible.');
end

% Eliminacion hacia adelante
for k = 1:n-1
    if A(k, k) == 0
        error('gauss_naive:pivoteNulo', 'Pivote nulo en k = %d. Use gauss_pivoteo.', k);
    end
    for i = k+1:n
        m = A(i, k) / A(k, k);
        A(i, k:n) = A(i, k:n) - m * A(k, k:n);
        b(i) = b(i) - m * b(k);
    end
end

x = sustitucion_atras(triu(A), b);
end
