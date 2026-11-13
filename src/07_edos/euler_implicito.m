function [t, Y] = euler_implicito(A, tspan, y0, n)
% EULER_IMPLICITO  Resuelve el sistema LINEAL y' = A*y con Euler implicito.
%
%   [t, Y] = euler_implicito(A, tspan, y0, n)
%
%   A     : matriz m x m constante (un escalar si hay una sola ecuacion)
%   tspan : [t0 tf]
%   y0    : condicion inicial (vector columna)
%   n     : numero de pasos
%
%   t : (n+1) x 1;  Y : (n+1) x m, una fila por instante.
%
%   Metodo (hacia atras):
%
%       y_{k+1} = y_k + h A y_{k+1}    =>    (I - h A) y_{k+1} = y_k
%
%   En cada paso se resuelve un sistema lineal. Como la matriz I - h A es la
%   misma en todos los pasos, se factoriza UNA sola vez (LU con pivoteo) y en
%   cada paso solo se hacen las sustituciones hacia adelante y hacia atras.
%   Es estable para sistemas rigidos (stiff).

if ~isnumeric(A)
    error('euler_implicito:matriz', 'A debe ser la matriz numerica del sistema y'' = A*y.');
end

y0 = y0(:);
m = numel(y0);
if ~isequal(size(A), [m, m])
    error('euler_implicito:dimension', 'A debe ser m x m, con m = numel(y0).');
end

t0 = tspan(1); tf = tspan(2);
h = (tf - t0) / n;
t = t0 + (0:n)' * h;
Y = zeros(n+1, m);
Y(1, :) = y0';

[L, U, P] = lu_factorizacion(eye(m) - h * A);    % P*(I - hA) = L*U

y = y0;
for k = 1:n
    y = sustitucion_atras(U, sustitucion_adelante(L, P * y));
    Y(k+1, :) = y';
end
end
