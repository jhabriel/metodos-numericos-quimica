function [t, Y] = heun(f, tspan, y0, n)
% HEUN  Resuelve y' = f(t,y) con el metodo de Heun (predictor-corrector).
%
%   [t, Y] = heun(f, tspan, y0, n)
%
%   Predictor (Euler):   y*_{k+1} = y_k + h f(t_k, y_k)
%   Corrector:           y_{k+1}  = y_k + h/2 [f(t_k,y_k) + f(t_{k+1}, y*_{k+1})]
%
%   Orden global 2. Mismos argumentos y salidas que euler_explicito.

t0 = tspan(1); tf = tspan(2);
h = (tf - t0) / n;
y0 = y0(:);
t = t0 + (0:n)' * h;
Y = zeros(n+1, numel(y0));
Y(1, :) = y0';
y = y0;
for k = 1:n
    f1 = f(t(k), y);
    yp = y + h * f1;                    % predictor
    f2 = f(t(k+1), yp);
    y = y + h / 2 * (f1 + f2);          % corrector
    Y(k+1, :) = y';
end
end
