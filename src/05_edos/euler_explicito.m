function [t, Y] = euler_explicito(f, tspan, y0, n)
% EULER_EXPLICITO  Resuelve y' = f(t,y) (escalar o sistema) con Euler explicito.
%
%   [t, Y] = euler_explicito(f, tspan, y0, n)
%
%   f     : handle f(t, y) -> vector columna
%   tspan : [t0 tf]
%   y0    : condicion inicial (vector columna)
%   n     : numero de pasos
%
%   t : (n+1) x 1;  Y : (n+1) x m, una fila por instante.
%
%   y_{k+1} = y_k + h f(t_k, y_k)

t0 = tspan(1); tf = tspan(2);
h = (tf - t0) / n;
y0 = y0(:);
t = t0 + (0:n)' * h;
Y = zeros(n+1, numel(y0));
Y(1, :) = y0';
y = y0;
for k = 1:n
    fy = f(t(k), y);
    y = y + h * fy(:);                  % fy(:) fuerza columna (por si f devuelve fila)
    Y(k+1, :) = y';
end
end
