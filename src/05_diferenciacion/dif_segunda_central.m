function d2 = dif_segunda_central(f, x, h)
% DIF_SEGUNDA_CENTRAL  Segunda derivada por diferencias centrales, O(h^2).
%
%   d2 = dif_segunda_central(f, x, h)
%
%   f''(x) ~ (f(x+h) - 2 f(x) + f(x-h)) / h^2

d2 = (f(x + h) - 2 * f(x) + f(x - h)) / h^2;
end
