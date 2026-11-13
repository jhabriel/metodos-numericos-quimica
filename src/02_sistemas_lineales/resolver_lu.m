function x = resolver_lu(A, b)
% RESOLVER_LU  Resuelve A*x = b usando la factorizacion LU propia.
%
%   x = resolver_lu(A, b)
%
%   P*A = L*U  =>  L*y = P*b (adelante),  U*x = y (atras).

[L, U, P] = lu_factorizacion(A);
y = sustitucion_adelante(L, P * b(:));
x = sustitucion_atras(U, y);
end
