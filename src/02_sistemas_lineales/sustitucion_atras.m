function x = sustitucion_atras(U, b)
% SUSTITUCION_ATRAS  Resuelve U*x = b con U triangular superior.
%
%   x = sustitucion_atras(U, b)

n = length(b);
b = b(:);
x = zeros(n, 1);
for i = n:-1:1
    if U(i, i) == 0
        error('sustitucion_atras:singular', 'Elemento diagonal nulo en la fila %d.', i);
    end
    s = 0;
    if i < n, s = U(i, i+1:n) * x(i+1:n); end
    x(i) = (b(i) - s) / U(i, i);
end
end
