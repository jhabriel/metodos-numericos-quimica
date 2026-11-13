function [a, r2, yajust] = regresion_lineal_general(x, y, base)
% REGRESION_LINEAL_GENERAL  Minimos cuadrados con modelo lineal en los parametros.
%
%   [a, r2, yajust] = regresion_lineal_general(x, y, base)
%
%   Modelo: y ~ a1*phi1(x) + a2*phi2(x) + ... + am*phim(x)
%
%   x, y : datos (vectores de igual longitud)
%   base : cell array de handles {@(x) ..., @(x) ...}, vectorizados
%
%   a      : coeficientes (m x 1), de las ecuaciones normales (Z'Z) a = Z'y
%   r2     : coeficiente de determinacion
%   yajust : valores ajustados en x
%
%   Ejemplos:
%     recta      : base = {@(x) ones(size(x)), @(x) x}
%     parabola   : base = {@(x) ones(size(x)), @(x) x, @(x) x.^2}
%     Arrhenius  : ln k = ln A - (Ea/R)*(1/T)  ->  x = 1/T, y = ln k, base de recta

x = x(:); y = y(:);
n = length(x);
m = numel(base);
if length(y) ~= n
    error('regresion_lineal_general:dimension', 'x e y deben tener la misma longitud.');
end
if n < m
    error('regresion_lineal_general:pocosDatos', 'Se necesitan al menos %d datos.', m);
end

Z = zeros(n, m);
for j = 1:m
    Z(:, j) = base{j}(x);
end

a = resolver_lu(Z' * Z, Z' * y);
yajust = Z * a;

St = sum((y - mean(y)).^2);
Sr = sum((y - yajust).^2);
if St == 0
    r2 = 1;
else
    r2 = 1 - Sr / St;
end
end
