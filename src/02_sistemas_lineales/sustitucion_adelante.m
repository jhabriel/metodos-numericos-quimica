function x = sustitucion_adelante(L, b)
% SUSTITUCION_ADELANTE  Resuelve L*x = b con L triangular inferior.
%
%   x = sustitucion_adelante(L, b)

n = length(b);
b = b(:);
x = zeros(n, 1);
for i = 1:n
    if L(i, i) == 0
        error('sustitucion_adelante:singular', 'Elemento diagonal nulo en la fila %d.', i);
    end
    s = 0;
    if i > 1, s = L(i, 1:i-1) * x(1:i-1); end
    x(i) = (b(i) - s) / L(i, i);
end
end
