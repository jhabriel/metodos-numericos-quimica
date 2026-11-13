function [t, Y] = euler_implicito(f, tspan, y0, n, J)
% EULER_IMPLICITO  Resuelve y' = f(t,y) con Euler implicito (hacia atras).
%
%   [t, Y] = euler_implicito(f, tspan, y0, n, J)
%
%   f     : handle f(t, y) -> vector columna
%   tspan : [t0 tf]
%   y0    : condicion inicial
%   n     : numero de pasos
%   J     : (opcional) handle del jacobiano df/dy, J(t, y). Si se omite
%           se aproxima por diferencias finitas.
%
%   y_{k+1} = y_k + h f(t_{k+1}, y_{k+1})
%   En cada paso se resuelve G(y) = y - y_k - h f(t_{k+1}, y) = 0 con
%   newton_sistema. Es estable para sistemas rigidos (stiff).

if nargin < 5, J = []; end

t0 = tspan(1); tf = tspan(2);
h = (tf - t0) / n;
y0 = y0(:);
m = numel(y0);
t = t0 + (0:n)' * h;
Y = zeros(n+1, m);
Y(1, :) = y0';
y = y0;
for k = 1:n
    tk1 = t(k+1);
    yk = y;
    G = @(z) z - yk - h * f(tk1, z);
    if isempty(J)
        JG = [];
    else
        JG = @(z) eye(m) - h * J(tk1, z);
    end
    y = newton_sistema(G, JG, yk + h * f(t(k), yk), 1e-12, 50);  % predictor Euler como inicio
    Y(k+1, :) = y';
end
end
