function [p, L] = lagrange(x, y)
% LAGRANGE  Coeficientes del polinomio interpolante de Lagrange.
%
%   [p, L] = lagrange(x, y)
%
%   ENTRADAS
%     x : abscisas de los nodos (n valores DISTINTOS)
%     y : valores en los nodos, y(i) = f(x(i)) (n valores)
%
%   SALIDAS
%     p : coeficientes del polinomio interpolante P(x), fila 1 x n, en
%         potencias DESCENDENTES (misma convencion que polyfit/polyval):
%
%             P(x) = p(1)*x^(n-1) + p(2)*x^(n-2) + ... + p(n-1)*x + p(n)
%
%     L : (opcional) matriz n x n con los coeficientes de los polinomios
%         base de Lagrange. La fila i contiene los coeficientes de L_i(x),
%         tambien en potencias descendentes, y se cumple p = y' * L.
%
%   METODO
%     El polinomio de grado <= n-1 que pasa por los n puntos (x_i, y_i) es
%
%         P(x) = sum_{i=1}^{n} y_i * L_i(x),
%
%         L_i(x) = prod_{j~=i} (x - x_j) / (x_i - x_j)
%
%     Cada L_i vale 1 en x_i y 0 en los demas nodos, por lo que
%     P(x_i) = y_i. El polinomio interpolante es unico. Aqui cada L_i se
%     expande multiplicando los factores lineales (x - x_j).
%
%   EJEMPLO
%     x = [0 1 2]; y = [1 3 7];        % puntos de P(x) = x^2 + x + 1
%     p = lagrange(x, y)               % devuelve [1 1 1]
%     polyval(p, 1.5)                  % evalua P(1.5) = 4.75
%
%     [p, L] = lagrange(x, y);
%     xx = linspace(0, 2, 100);
%     plot(xx, polyval(L(1, :), xx))   % grafica la base L_1(x)
%
%   NOTAS
%     - Para evaluar usar polyval(p, xq).
%     - Interpolar es fiable dentro del rango de los nodos; fuera de el
%       (extrapolacion) o con muchos nodos equiespaciados el error puede
%       crecer mucho (fenomeno de Runge). Un grado alto no implica mejor
%       aproximacion.
%     - Con muchos nodos o abscisas grandes, los coeficientes en potencias
%       de x se vuelven sensibles a errores de redondeo (el problema esta
%       mal condicionado). Conviene centrar/escalar x antes de interpolar.

x = x(:); y = y(:);
n = length(x);
if length(y) ~= n
    error('lagrange:dimension', 'x e y deben tener la misma longitud.');
end
if numel(unique(x)) ~= n
    error('lagrange:nodosRepetidos', 'Los nodos x deben ser distintos.');
end

L = zeros(n, n);
for i = 1:n
    Li = 1;                  % polinomio constante 1
    den = 1;                 % producto de (x_i - x_j)
    for j = [1:i-1, i+1:n]
        Li = conv(Li, [1, -x(j)]);      % multiplica por (x - x_j)
        den = den * (x(i) - x(j));
    end
    L(i, :) = Li / den;      % L_i(x_i) = 1
end

p = y' * L;                  % P = sum_i y_i * L_i
end
