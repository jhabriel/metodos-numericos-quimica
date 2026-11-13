function [xmin, fmin, iter] = seccion_aurea(f, a, b, tol, maxit)
% SECCION_AUREA  MINIMO de f en [a,b] (unimodal) por seccion dorada.
%
%   [xmin, fmin, iter] = seccion_aurea(f, a, b, tol, maxit)
%
%   ATENCION: esta funcion solo MINIMIZA. Para hallar un MAXIMO de f hay que
%   minimizar g(x) = -f(x) y luego cambiar el signo del resultado:
%
%       [xmax, gmin] = seccion_aurea(@(x) -f(x), a, b);
%       fmax = -gmin;          % el maximo de f es el opuesto del minimo de g
%
%   ENTRADAS
%     f     : handle de la funcion, unimodal en [a,b] (un solo minimo)
%     a, b  : extremos del intervalo
%     tol   : tolerancia sobre el ancho del intervalo (defecto 1e-6)
%     maxit : maximo de iteraciones (defecto 200)
%
%   SALIDAS
%     xmin  : abscisa del minimo
%     fmin  : valor f(xmin)
%     iter  : iteraciones realizadas
%
%   METODO
%     Se mantienen dos puntos interiores x1 < x2 colocados con la razon
%     aurea R = (sqrt(5)-1)/2. Se descarta el subintervalo que no contiene
%     al minimo y se reutiliza uno de los puntos, por lo que cada iteracion
%     solo evalua f una vez.

if nargin < 4, tol = 1e-6; end
if nargin < 5, maxit = 200; end

R = (sqrt(5) - 1) / 2;          % 0.6180...
x1 = b - R * (b - a);
x2 = a + R * (b - a);
f1 = f(x1); f2 = f(x2);

for iter = 1:maxit
    if (b - a) < tol, break; end
    if f1 < f2                  % el minimo esta en [a, x2]
        b = x2;
        x2 = x1; f2 = f1;
        x1 = b - R * (b - a);
        f1 = f(x1);
    else                        % el minimo esta en [x1, b]
        a = x1;
        x1 = x2; f1 = f2;
        x2 = a + R * (b - a);
        f2 = f(x2);
    end
end

xmin = (a + b) / 2;
fmin = f(xmin);
end
